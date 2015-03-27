% MEGSearchlightRDMs_source
%
% Cai Wingfield 2015-03

function [RDMsPaths] = MEGSearchlightRDMs_source(meshPaths, slMasks, adjacencyMatrix, STCMetadata, userOptions)

import rsa.*
import rsa.meg.*
import rsa.util.*

%% Constants

usingMasks = ~isempty(userOptions.maskNames);

nSubjects = numel(userOptions.subjectNames);

%% File paths

returnHere = pwd; % We'll come back here later

RDMsDir = fullfile(userOptions.rootPath, 'RDMs');

file_i = 1;
for subject_i = 1:nSubjects
    thisSubjectName = userOptions.subjectNames{subject_i};
    for chi = 'LR'
        if usingMasks
            RDMsFile = ['searchlightRDMs_masked_', thisSubjectName, '-' lower(chi) 'h'];
        else
            RDMsFile = ['searchlightRDMs_',        thisSubjectName, '-' lower(chi) 'h'];
        end
        RDMsPaths(subject_i).(chi) = fullfile(RDMsDir, RDMsFile);
        
        % We'll check all the files to be saved to see if they have already
        % been saved.
        promptOptions.checkFiles(file_i).address = RDMsPaths(subject_i).(chi);
        file_i = file_i + 1;
    end
end

promptOptions.functionCaller = 'MEGSearchlightRDMs_source';
promptOptions.defaultResponse = 'S';

overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
%% Apply searchlight
    
parfor subject_i = 1:nSubjects
    thisSubjectName = userOptions.subjectNames{subject_i};

    % Work on each hemisphere separately
    for chi = 'LR'
        
        % We'll only do the searchlight if we haven't already done so,
        % unless we're being told to overwrite.
        if exist(RDMsPaths(subject_i).(chi), 'file') && ~overwriteFlag
            prints('Searchlight already performed in %sh hemisphere of subject %d. Skipping.', lower(chi), subject_i);
        else
            prints('Shining RSA searchlight in the %sh source mesh of subject %d of %d (%s)...', lower(chi), subject_i, nSubjects, thisSubjectName);
            
            single_hemisphere_searchlight( ...
                STCMetadata, ...
                meshPaths(subject_i).(chi), ...
                RDMsPaths(subject_i).(chi), ...
                RDMsDir, ...
                slMasks([slMasks.chi] == chi), ...
                adjacencyMatrix, ...
                userOptions);

            %% Done
            prints('Done with subeject %d''s %sh side.', subject_i, lower(chi));

        end%if:overwrite
    end%for:chi
end%for:subject

cd(returnHere); % And go back to where you started

end%function

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

% Computes and saves searchlight RDMs for a single hemisphere of a single
% subject.
function single_hemisphere_searchlight(STCMetadata, meshPath, RDMsPath, RDMsDir, slMask, adjacencyMatrix, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.stat.*
    import rsa.util.*

    [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions);

    maskedMeshes = directLoad(meshPath, 'sourceMeshes');

    searchlightRDMs = searchlightMappingRDMs_MEG_source(maskedMeshes, slMask, adjacencyMatrix, slSpec, userOptions); %#ok<NASGU>

    %% Saving RDM maps

    prints('Saving data RDMs to %s.', RDMsPath);
    gotoDir(RDMsDir);
    save('-v7.3', RDMsPath, 'searchlightRDMs');
end%function

% Computes and returns searchlight RDMs for a single hemisphere for a
% single subject.
%
% Based on Li Su's script
% Cai Wingfield 2015-03
function [searchlightRDMs] = searchlightMappingRDMs_MEG_source(maskedMeshes, indexMask, adjacencyMatrix, slSpec, userOptions)

import rsa.*
import rsa.meg.*
import rsa.rdm.*
import rsa.stat.*
import rsa.util.*

[nVertices_masked, nTimePoints_data, nConditions, nSessions] = size(maskedMeshes);

% The number of positions the sliding window will take.
nWindowPositions = size(slSpec.windowPositions, 1);

%% map the volume with the searchlight

% Preallocate looped matrices for speed
searchlightRDMs(1:nVertices_masked, 1:nWindowPositions) = struct('RDM', size(vectorizeRDM(zeros(nConditions))));

% For display purposes only
nVertsSearched = 0;

% Search the vertices
for v_i = 1:numel(indexMask.vertices)
    % v_i loops through the *indices of* vertices in the mask
    % v is the vertex itself
    
    v = indexMask.vertices(v_i);
    
    % Determine which vertexes are within the radius of the currently-picked vertex
    searchlightPatchVs = [v, adjacencyMatrix(v,:)];
    
    % Restrict to verticies inside mask.
    % This also removes any nans.
    searchlightPatchVs = intersect(searchlightPatchVs, indexMask.vertices);
    
    % Now we need to convrt the vertices into vertex *indices*.  For
    % example, our mask may be vertices [5, 6, 7], but since there will
    % only be three datapoints inside each of the masked meshes, we need to
    % be able to refer to vertex 1, 2 and 3.
    
    % Unfortunately, the only way I can think to do this is with a reverse
    % lookup loop, which isn't too efficient.  Hopefuly the small
    % searchlight patch will be fast enough, though.
    searchlightPatchVIs = [];
    for slVertex = searchlightPatchVs
        slVertex_i = find(searchlightPatchVs == slVertex);
        searchlightPatchVIs = [searchlightPatchVIs, slVertex_i];
    end%for
    
    % Search through time
    window_i = 0;
    for window = slSpec.windowPositions'
        % thisWindow is the indices of timepoints in each window
        thisWindow = window(1):window(2);
        window_i = window_i + 1;
        
        searchlightPatchData = maskedMeshes(searchlightPatchVIs, thisWindow, :, :); % (vertices, time, condition, session)
        
        if ~userOptions.regularized
            
            switch lower(userOptions.searchlightPatterns)
                case 'spatial'
                    % Spatial patterns: median over time window
                    searchlightPatchData = median(searchlightPatchData, 2); % (vertices, 1, conditions, sessions)
                    searchlightPatchData = squeeze(searchlightPatchData); % (vertices, conditions, sessions);
                case 'temporal'
                    % Temporal patterns: mean over vertices within searchlight
                    searchlightPatchData = mean(searchlightPatchData, 1); % (1, timePoints, conditions, sessions)
                    searchlightPatchData = squeeze(searchlightPatchData); % (timePionts, conditions, sessions)
                case 'spatiotemporal'
                    % Spatiotemporal patterns: all the data concatenated
                    searchlightPatchData = reshape(searchlightPatchData, [], size(searchlightPatchData, 3), size(searchlightPatchData, 4)); % (dataPoints, conditions, sessions)
            end%switch:userOptions.sensorSearchlightPatterns
            
            % Average RDMs over sessions
            searchlightRDM_square = zeros(nConditions);
            for session = 1:nSessions
                sessionRDM = squareform(pdist(squeeze(searchlightPatchData(:,:,session))',userOptions.distance));
                searchlightRDM_square = searchlightRDM_square + sessionRDM;
            end%for:sessions
            searchlightRDM_square = searchlightRDM_square / nSessions;
            
        else
            % TODO: Look into this closer
            % data regularization based on algorithm by Diedrichson et al 2011 - updated 12-12 IZ
            tempMesh = reshape(searchlightPatchData, [], size(searchlightPatchData, 3), size(searchlightPatchData, 4));
            searchlightPatchData = zeros(size(tempMesh, 1), size(tempMesh, 2) * size(tempMesh, 3)); % (data, conditions, sessions)
            
            % combining session-wise trials
            kk = 1;
            for j = 1:size(tempMesh,2)
                for i = 1:nSessions
                    searchlightPatchData(:, kk) = (tempMesh(:, j, i));
                    kk = kk + 1;
                end
            end
            
            r_matrix = g_matrix(zscore(squeeze(searchlightPatchData(:,:)))', nConditions, size(currentTimeWindow,2));
            searchlightRDM_square = (1 - r_matrix);
            
            if isnan(searchlightRDM_square) % sessions and conditions should be optimal
                error('Cannot calculate covariance matrix. Try reducing number of conditions');
            end
        end
        
        % Store results to be retured.
        searchlightRDMs(v_i, window_i).RDM = vectorizeRDM(searchlightRDM_square);
        
    end%for:window
    
    % Indicate progress every once in a while...
    nVertsSearched = nVertsSearched + 1;
    if mod(nVertsSearched, 100) == 0
        prints('%d vertices searched', nVertsSearched);
    end%if
    
end%for:v

end%function

