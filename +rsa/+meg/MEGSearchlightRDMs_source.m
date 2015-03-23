% MEGSearchlightRDMs_source
%
% Cai Wingfield 2015-03

function [RDMsPath] = MEGSearchlightRDMs_source(subject_i, chi, maskedMeshes, slMask, adjacencyMatrix, STCMetadata, userOptions)

import rsa.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

subjectName = userOptions.subjectNames{subject_i};

usingMasks = ~isempty(userOptions.maskNames);

%% File paths

RDMsDir = fullfile(userOptions.rootPath, 'RDMs');

if usingMasks
    RDMsFile = ['searchlightRDMs_masked_', subjectName, '-' lower(chi) 'h'];
else
    RDMsFile = ['searchlightRDMs_',        subjectName, '-' lower(chi) 'h'];
end
RDMsPath = fullfile(RDMsDir, RDMsFile);

promptOptions.functionCaller = 'MEGSearchlightRDMs_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = RDMsPath;

overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
%% Apply searchlight

if overwriteFlag
    
    prints('Shining RSA searchlight in the source mesh of subject %d of %d...', subject_i, numel(userOptions.subjectNames));
    
    tic;%1
    
    [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions);
    
    searchlightRDMs = rsa.meg.searchlightMappingRDMs_MEG_source(maskedMeshes, slMask, adjacencyMatrix, slSpec, userOptions); %#ok<NASGU>

    %% Saving RDM maps
    
    prints('Saving data RDMs to %s.', RDMsPath);
    gotoDir(RDMsDir);
    save('-v7.3', RDMsPath, 'searchlightRDMs');
    
    %% Done
    t = toc;%1
    prints('That took %s seconds.', t);
else
    prints('Searchlight already applied for subject %s, %s side, skipping it.', subjectName, lower(chi));
end

cd(returnHere); % And go back to where you started

end%function

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

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
    if mod(nVertsSearched, 500) == 0
        prints('%d vertices searched', nVertsSearched);
    end%if
    
end%for:v

end%function

