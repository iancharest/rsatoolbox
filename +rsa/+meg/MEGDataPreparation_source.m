% MEGDataPreparation_source is a function designed to take MEG data and load it
% into the correct format for the rest of the toolbox to use.
%
% It is based on Su Li's code
%
% [meshPaths, STCMetadata] = MEGDataPreparation_source(
%                                          betas,
%                                          userOptions,
%                                         ['mask', indexMask,]
%                                         ['subject_i', subject_i,]
%                                         ['chi', 'L'|'R'])
%
%       betaCorrespondence --- The array of beta filenames.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image.
%
%       userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.betaPath
%                        A string which contains the absolute path to the
%                        location of the beta images. It can contain the
%                        following wildcards which would be replaced as
%                        indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%                                [[betaIdentifier]]
%                                        To be replaced by filenames as provided
%                                        by betaCorrespondence.
%
%       indexMask --- a struct containing at least:
%                indexMask.vertices
%                        A vector of integers which represent the vertices
%                        inside the mask.
%
%       subject_i --- a subject index
%                If not provided, all subjects will be looped through.
%
%       chi --- 'L' or 'R'
%                If not provided, both hemispheres will be looped through.
%
%       sourceMeshes
%                sourceMeshes.(subjectName).(L|R)(vertices, timepoints, condition, session)
%                This data is downsampled, as per userOptions preferences.
%
% REQUIRES MATLAB VERSION 7.3 OR LATER
%
% Cai Wingfield 2010-05, 2010-06, 2015-03
% updated by Li Su 2-2012
% updated by Fawad 3-12014

function [meshPaths, STCMetadata] = MEGDataPreparation_source(betas, userOptions, varargin)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% Parse inputs

% 'masks'
nameMask = 'mask';
% We don't do any validation checks on the mask
checkMask = @(x) (true);
defaultMask = {};

% 'subject_i'
nameSubjectI    = 'subject_i';
checkSubjectI   = @(x)(isnumeric(x) && x >= 0 && x <= numel(userOptions.subjectNames));
defaultSubjectI = 0;

% 'chi'
nameChi = 'chi';
validChis = {'L', 'R', ''};
checkChi = @(x) (any(validatestring(x, validChis)));
defaultChi = '';

% Set up parser
ip = inputParser;
ip.CaseSensitive = false;
ip.StructExpand = false;

% Parameters
addParameter(ip, nameMask, defaultMask, checkMask);
addParameter(ip, nameSubjectI, defaultSubjectI, checkSubjectI);
addParameter(ip, nameChi, defaultChi, checkChi);

% Parse the inputs
parse(ip, varargin{:});

% If subject_i was not given a default value, then it will be a positive
% integer, and that's the subject_i we'll use
singleSubject = (ip.Results.subject_i > 0);

% If chi was not given a default value, then it will be either 'L' or 'R',
% and we'll use that.
singleHemisphere = ~isempty(ip.Results.chi);

% If masks was given a default value, we're not going to do any masking.
usingMask = ~isempty(ip.Results.mask);

% We'll return to the pwd when the function has finished
returnHere = pwd;

% If we're running a single subject, the we "loop" once on that
% subject. If we're running all subjects, we loop through them all
if singleSubject
    firstSubject_i = ip.Results.subject_i;
    lastSubject_i = ip.Results.subject_i;
else
    firstSubject_i = 1;
    lastSubject_i = numel(userOptions.subjectNames);
end

% If we're running a single hemisphere, then we "loop" once on that
% side. If we're running both sides, we loop through them both.
if singleHemisphere
    chis = ip.Results.chi;
else
    chis = 'LR';
end

[nSessions, nConditions] = size(betas);

% Before we've read any data, we don't know the following values, which
% will be set appropriately when we load the first piece of data.
% It's probably not necessary to do this, but I can't deal with Matlab's
% gross lack of scoping, so by god I'll do it anyway.
dataEverRead            = false;
nVertices_raw           = NaN;
nTimepoints_raw         = NaN;
nTimepoints_downsampled = NaN;
STCMetadata             = struct();

% Some conditions will have been rejected, and we'll record those in
% this text file.
imageDataPath = fullfile(userOptions.rootPath, 'ImageData');
missingFilesLog = fullfile(imageDataPath, 'missingFilesLog.txt');

% Where data will be saved
gotoDir(imageDataPath);

%% Loop over all subjects under consideration
for subject_i = firstSubject_i:lastSubject_i
    
    % Figure out the subject's name
    thisSubjectName = userOptions.subjectNames{subject_i};
    
    %% Loop over all hemispheres under consideration
    for chi = chis
        
        % Loop over sessions and conditions
        for session_i = 1:nSessions
            for condition_i = 1:nConditions
                
                % Then read the brain data (for this session, condition)
                readPath = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', betas(session_i, condition_i).identifier, '[[subjectName]]', thisSubjectName, '[[LR]]', lower(chi));
                
                dataReadSuccessfully = false;
                try
                    MEGData_stc = mne_read_stc_file1(readPath);
                    dataReadSuccessfully = true;
                catch ex
                    % when a trial is rejected due to artifact, this item
                    % is replaced by NaNs. Li Su 3-2012
                    prints(['Warning: Failed to read data for condition ' num2str(condition_i) '... Writing NaNs instead.']);
                    % Log the missing file
                    dlmwrite(missingFilesLog, str2mat(replaceWildcards(betas(session_i, condition_i).identifier, '[[subjectName]]', thisSubjectName)), 'delimiter', '', '-append');
                end
                
                if dataReadSuccessfully
                    
                    %% First time data is read
                    
                    % If this was the first time we read data, we can
                    % record the sizes and do the preallocation
                    if ~dataEverRead
                        
                        % We only want to do this once so we'll
                        % remember now that we've read the data and
                        % don't ned to do this again.
                        dataEverRead = true;
                        
                        % Raw data sizes, before downsampling
                        
                        % The number of vertices and timepoints in the
                        % raw data
                        [nVertices_raw, nTimepoints_raw] = size(MEGData_stc.data);
                        % The time index of the first datapoint, in
                        % seconds.
                        firstDatapointTime_raw = MEGData_stc.tmin;
                        % The interval between successive datapoints in
                        % the raw data, in seconds.
                        timeStep_raw = MEGData_stc.tstep;
                        
                        %% Downsampling constants
                        
                        % Sanity checks for downsampling
                        if nVertices_raw < userOptions.targetResolution
                            error('MEGDataPreparation_source:InsufficientSpatialResolution', 'There aren''t enough vertices in the raw data to meet the target resolution.');
                        end
                        
                        % Now the actual downsampling targets can be
                        % calculated
                        
                        % We will be downsampling in space and time, so
                        % we calculate some useful things here.
                        % The number of timepoints in the downsampled
                        % data
                        nTimepoints_downsampled = numel(1:userOptions.temporalDownsampleRate:nTimepoints_raw);
                        timeStep_downsampled = timeStep_raw * userOptions.temporalDownsampleRate;
                        
                        % Time time index of the first datapoint
                        % doesn't change in the downsampled data
                        firstDatapointTime_downsampled = firstDatapointTime_raw;
                        
                        % The time index of the last datapoint may
                        % change in the downsampled data, so should be
                        % recalculated
                        lastDatapointTime_downsampled = firstDatapointTime_downsampled + (nTimepoints_downsampled * timeStep_downsampled);
                        
                        % This metadata struct will be useful for
                        % writing appropriate files in future. This new
                        % metadata should reflect the resolution and
                        % specifices of the data which
                        % MEGDataPreparation_source produces.
                        STCMetadata.tmin     = firstDatapointTime_downsampled;
                        STCMetadata.tmax     = lastDatapointTime_downsampled;
                        STCMetadata.tstep    = timeStep_downsampled;
                        
                        %% Apply mask
                        
                        % If we're using masks, we only want to include
                        % those vertices which are inside the mask.
                        if usingMask
                            STCMetadata.vertices = sort(ip.Results.mask.vertices);
                        else
                            % If we're not using a mask, we still need to
                            % downsample the mesh to the target resolution.
                            % Luckily, downsampling is just a matter of
                            % taking low-numbered vertices, due to the way
                            % they are laid out.
                            STCMetadata.vertices = 1:userOptions.targetResolution;
                        end%if
                    end
                    
                    %% Every time data is read
                    
                    % Store the data in the mesh, masking and downsampling
                    % as we go.
                    sourceMeshes(:, :, condition_i, session_i) = ...
                        MEGData_stc.data( ...
                            ...% Downsample and mask space
                            STCMetadata.vertices, ...
                            ...% Downsample time by subsampling the timepoints
                            1:userOptions.temporalDownsampleRate:end); % (vertices, time, condition, session)
                else
                    % Make sure it actually has NaNs in if there was an
                    % error for this condition
                    sourceMeshes(:, :, condition_i, session_i) = NaN(numel(STCMetadata.vertices), nTimepoints_downsampled);
                end
            end%for:condition
        end%for:session
        
        % TODO: This won't work if more than one subject or chi are looped
        % TODO: over.  Do we want to prevent this entirely? Or change
        % TODO: what's returned depending on the optional arguments?
        meshPaths = fullfile(imageDataPath, [userOptions.analysisName, '_', thisSubjectName, '_CorticalMeshes.mat']);
        
        gotoDir(imageDataPath);
        save('-v7.3', meshPaths, 'sourceMeshes');
        
        prints('Subject %d''s %s-hemisphere data read successfully!', subject_i, chi);
        dlmwrite(missingFilesLog, '', '-append');
    end%for:chi
end%for:subjects

cd(returnHere); % Go back
