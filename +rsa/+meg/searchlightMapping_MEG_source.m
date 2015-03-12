% [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_source(singleMesh, model, mask, userOptions, localOptions)
%
% Based on Li Su's script
% CW 2010-05, 2015-03
% updated by Li Su 3-2012

function [smm_rs, searchlightRDMs] = searchlightMapping_MEG_source(singleSubjectMesh, indexMask, modelRDM, partialModelRDMs, adjacencyMatrix, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% Get parameters

modelRDM_utv = squeeze(unwrapRDMs(vectorizeRDMs(modelRDM)));

if userOptions.partial_correlation
    % TODO: should be transposed?
    control_for_modelRDMs = unwrapRDMs(vectoriseRDMs(partialModelRDMs));
end

[nVertices, epochLength, nConditions, nSessions] = size(singleSubjectMesh);

% Number of DATA points to loop with given the width and time step of
% searchlight updated by IZ 09-12
nTimePoints = floor((epochLength - (userOptions.temporalSearchlightWidth * userOptions.toDataPoints)) / ...
    (userOptions.temporalSearchlightResolution * userOptions.toDataPoints * userOptions.temporalDownsampleRate));


%% similarity-graph-map the volume with the searchlight

% Preallocate looped matrices for speed
smm_rs = zeros([nVertices, nTimePoints, nModels]);

for v = indexMask.vertices
    
    % Determine which vertexes are within the radius of the currently-picked vertex
    verticesCurrentlyWithinRadius = [v, adjacencyMatrix(v,:)];
    
    % Restrict to verticies inside mask.
    % This also removes any nans.
    % All searchlight run as masks, including full-brain searchlights (update IZ 03/12)
    verticesCurrentlyWithinRadius = intersect(verticesCurrentlyWithinRadius, vertices);
    
    % TODO: Why +1 ?
    for t = 1:nTimePoints+1
        
        % Work out the current time window
        % converted to data points - updated by IZ 09-12
        currentTimeStart = (t - 1) * ...
            (userOptions.temporalSearchlightResolution * userOptions.temporalDownsampleRate * userOptions.toDataPoints) + 1;
        currentTimeWindow = ceil((currentTimeStart : currentTimeStart + ...
            (userOptions.temporalSearchlightWidth * userOptions.toDataPoints) - 1));
        
        searchlightPatchData = singleSubjectMesh(verticesCurrentlyWithinRadius, currentTimeWindow, :, :); % (vertices, time, condition, session)
        
        % Average across sessions
        
        if not(userOptions.regularized)
            
            % Median over the time window
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
            
            % Preallocate
            searchlightRDM = zeros(nConditions);
            
            for session = 1:userOptions.nSessions
                searchlightRDM = searchlightRDM + squareform(pdist(squeeze(searchlightPatchData(:,:,session))',userOptions.distance));
            end%for:sessions
            
        else
            % data regularization based on algorithm by Diedrichson et al 2011 - updated 12-12 IZ
            tempMesh = reshape(searchlightPatchData, [], size(searchlightPatchData, 3), size(searchlightPatchData, 4));
            searchlightPatchData = zeros(size(tempMesh, 1), size(tempMesh, 2) * size(tempMesh, 3)); % (data, conditions, sessions)
            
            % combining session-wise trials
            kk = 1;
            for j = 1:size(tempMesh,2)
                for i = 1:userOptions.nSessions
                    searchlightPatchData(:, kk) = (tempMesh(:, j, i));
                    kk = kk + 1;
                end
            end
            
            r_matrix = g_matrix(zscore(squeeze(searchlightPatchData(:,:)))', userOptions.nConditions, size(currentTimeWindow,2));
            searchlightRDM = searchlightRDM + (1 - r_matrix);
            
            if isnan(searchlightRDM) % sessions and conditions should be optimal
                error('Cannot calculate covariance matrix. Try reducing number of conditions');
            end
        end
        
        searchlightRDM = searchlightRDM / userOptions.nSessions;
        
        searchlightRDM = vectorizeRDM(searchlightRDM);
        
        % Locally store the full brain's worth of indexed RDMs.
        searchlightRDMs(v, t).RDM = searchlightRDM;
        
        % TODO: Refactor this into general method so it can be used
        % TODO: anywhere
        if strcmpi(userOptions.distanceMeasure, 'Kendall_taua')
            rs = rankCorr_Kendall_taua(searchlightRDM', modelRDM_utv');
        elseif userOptions.partial_correlation
            % TODO: Consider partialcorr with kendall's tau
            rs = partialcorr(searchlightRDM', modelRDM_utv', control_for_modelRDMs', 'type', userOptions.distanceMeasure, 'rows','pairwise');
        else
            rs = corr(searchlightRDM', modelRDM_utv', 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
        end
        
        smm_rs(v, t, :) = rs;
        
    end%for:t
    
    % Indicate progress every once in a while...
    if mod(v, floor(length(vertices) / 20)) == 0, fprintf('.'); end%if
    
end%for:vertices

if userOptions.fisher
    smm_rs = fisherTransform(smm_rs);
end%if

end%function
