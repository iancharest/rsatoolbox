% MEGDataMasking_source is a module which masks source-localised EMEG
% time-course data to RoIs specified by lists of time- and vertex-indices.
%
% maskedMeshes = MEGDataMasking_source(
%                                      sourceMeshes,
%                                      indexMasks,
%                                      betaCorrespondence,
%                                      userOptions
%                                     )
%
%        sourceMeshes --- The unmasked sourceMeshes topographies.
%                Contains the subject source-reconstructed mesh data in a
%                by-hemisphere struct:
%                sourceMeshes.L(vertices, timepoints, condition, session)
%                sourceMeshes.R(vertices, timepoints, condition, session)
%
%        indexMasks --- The specification of the masks which will be applied
%                       to the data.
%                indexMasks.vertices
%                        A vector of indices inside the mask (indices above
%                       10242 are ignored).
%                indexMasks.chi
%                        Either "L" or "R", depending on which hemisphere the
%                       mask indices refer to. Cross-hemisphere RoIs aren't
%                       currently supported...
%                indexMasks.name
%                        The name of the mask (hyphens replaced by
%                        underscores.
%
%        betaCorrespondence --- The array of beta filenames.
%                betas(condition, session).identifier is a string which referrs
%                to the filename (not including path) of the SPM beta image.
%                (Or, if not using SPM, just something, as it's used to
%                determine the number of conditions and sessions.)
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in sensorImages.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of indexMasks.
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_MaskedSensors.mat
%                        Contains the raw timecourses in a structure such that
%                        maskedSensors.masked.(subjectName) is a [nChannels
%                        nConditions nSessions]-sized matrix (median across time
%                        window).
%        userOptions.rootPath/Details/
%                userOptions.analysisName_MEGDataMasking_sensor_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%  
% Cai Wingfield 2010-09, 2015-03
% Updated Isma Zulfiqar 11-2012
% Updated Fawad 26/03/2014

function [maskedMeshes] = MEGDataMasking_source(sourceMeshes, indexMasks, betaCorrespondence, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

maximumVertexIndex = userOptions.targetResolution;

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('MEGDataMasking_source:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('MEGDataMasking_source:NoRootPath', 'rootPath must be set. See help'); end%if

% The analysisName will be used to label the files which are eventually saved.
MaskedBrainsFilename = [userOptions.analysisName, '_MaskedMeshes.mat'];
DetailsFilename = [userOptions.analysisName, '_MEGDataMasking_source_Details.mat'];

promptOptions.functionCaller = 'MEGDataMasking_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag

	% Data
	nMasks = numel(userOptions.maskNames);

	%% Get Data

	nSubjects = numel(userOptions.subjectNames);
    [nVertices, nTimepoints, nConditions, nSessions] = size(sourceMeshes.L);

	for mask_i = 1:nMasks
	
		% Which mask is this?
        thisMaskName = indexMasks(mask_i).name;
        
        % Load the mask data into a vector
        maskIndices = indexMasks(mask_i).vertices;
        maskIndices = maskIndices(maskIndices <= maximumVertexIndex);
        % TODO: This is confusing
        timeIndices = userOptions.maskTimetoDataPoints.(thisMaskName); % updated IZ 11-12
        % timeIndices = indexMasks.(thisMask).(thisTimeWindow).timeIndices;
        chi = indexMasks(mask_i).chi;

        for subject_i = 1:nSubjects % and for each subject...

            % Figure out which subject this is
            thisSubjectName = userOptions.subjectNames{subject_i};

            % Get the mesh for this subject
            thisMesh = sourceMeshes(subject_i).(chi);

            % Mask the data
            if nSessions == 1
                maskedMesh = thisMesh(maskIndices, timeIndices(1):timeIndices(2), :); % (vertices, timePointes, conditions) % updated IZ 11-12
            else
                maskedMesh = thisMesh(maskIndices, timeIndices(1):timeIndices(2), :, :); % (vertices, timePointes, conditions, sessions) % updated IZ 11-12
            end
         if ~userOptions.regularized
            % Reduce to correct data type
            switch lower(userOptions.searchlightPatterns)
                % For spatial patterns, median across the time window
                case 'spatial'
                    reducedMaskedMesh = zeros(size(maskedMesh, 1), size(maskedMesh, 3), size(maskedMesh, 4)); % (data, conditions, sessions)
                    reducedMaskedMesh(:,:,:) = median(maskedMesh, 2);
                % For temporal patterns, average across the vertices.
                case 'temporal'
                    reducedMaskedMesh = zeros(size(maskedMesh, 2), size(maskedMesh, 3), size(maskedMesh, 4)); % (data, conditions, sessions)
                    reducedMaskedMesh(:,:,:) = mean(maskedMesh, 1);
                % For spatiotemporal patterns, concatenate the vertex timecourses
                case 'spatiotemporal'
                    reducedMaskedMesh = reshape(maskedMesh, [], size(maskedMesh, 3), size(maskedMesh, 4)); % update IZ 11-12
            end % switch

          else% update IZ 11-12 
             % case 'regularized' 
             tempMesh = reshape(maskedMesh, [], size(maskedMesh, 3), size(maskedMesh, 4));   
             reducedMaskedMesh = zeros(size(maskedMesh, 1)*size(maskedMesh,2), size(maskedMesh, 3)* size(maskedMesh, 4)); % (data, conditions, sessions)

              % combining session-wise trials
             k=1;
             for j=1:size(tempMesh,2)
                for i=1:nSessions
                    reducedMaskedMesh(:,k) = (tempMesh(:,j,i)); 
                    k=k+1;
                end
             end
          end % if regularized                             

            % Store in struct
            maskedMeshes.(thisMaskName).(thisSubjectName) = reducedMaskedMesh;
            % maskedMeshes.(thisMask) = reducedMaskedMesh; % IZ 11/12

        end%for:subject
	end%for:mask

	%% Save relevant info

	timeStamp = datestr(now);
    
	fprintf(['Saving masked data to ' fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'ImageData');
	save(MaskedBrainsFilename, 'maskedMeshes');

	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');%, 'searchlightPatterns');
	
else
	fprintf(['Loading previously saved RoIs from ' fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename) '...\n']);
	maskedMeshes = load(fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename), 'maskedMeshes');
end%if

cd(returnHere); % And go back to where you started
