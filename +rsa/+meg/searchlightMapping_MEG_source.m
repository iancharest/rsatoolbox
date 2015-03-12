% [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_source(singleMesh, model, mask, userOptions, localOptions)
%
% Based on Li Su's script
% CW 2010-05, 2015-03
% updated by Li Su 3-2012

function [smm_rs, smm_ps, searchlightRDMs] = searchlightMapping_MEG_source(singleSubjectMesh, modelRDM, partialModelRMDs, adjacencyMatrix, userOptions)

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

% vertices change on the basis of maskng flag's value IZ 11-12
% updated: all searchlight run as masks IZ 03/12
vertices = userOptions.maskIndices.(userOptions.chi);

% Preallocate looped matrices for speed
smm_ps = zeros([nVertices, nTimePoints, nModels]);
smm_rs = zeros([nVertices, nTimePoints, nModels]);

for vertex = vertices
    
    % Determine which vertexes are within the radius of the currently-picked vertex
    verticesCurrentlyWithinRadius = adjacencyMatrix(vertex,:);
    
    % remove nans
    verticesCurrentlyWithinRadius = verticesCurrentlyWithinRadius(~isnan(verticesCurrentlyWithinRadius));
    
    % add current vertex
    verticesCurrentlyWithinRadius = [vertex, verticesCurrentlyWithinRadius]; % add by Li Su 1-02-2010
    
    % If masks are used, finding corresponding mask indices - update IZ 11-12
    if userOptions.maskingFlag
        verticesCurrentlyWithinRadius = intersect(verticesCurrentlyWithinRadius, vertices);
    end
    
    for t = 1:nTimePoints+1
        
        % Work out the current time window
        % converted to data points - updated by IZ 09-12
        currentTimeStart = (t - 1) * ...
            (userOptions.temporalSearchlightResolution * userOptions.temporalDownsampleRate * userOptions.toDataPoints) + 1;
        currentTimeWindow = ceil((currentTimeStart : currentTimeStart + ...
            (userOptions.temporalSearchlightWidth * userOptions.toDataPoints) - 1));
        
        currentData = singleSubjectMesh(verticesCurrentlyWithinRadius, currentTimeWindow, :, :); % (vertices, time, condition, session)
        
        searchlightRDM = zeros(nConditions, nConditions);
        
        % Average across sessions
        
        if not(userOptions.regularized)
            
            % Median over the time window
            switch lower(userOptions.searchlightPatterns)
                case 'spatial'
                    % Spatial patterns: median over time window
                    currentData = median(currentData, 2); % (vertices, 1, conditions, sessions)
                    currentData = squeeze(currentData); % (vertices, conditions, sessions);
                case 'temporal'
                    % Temporal patterns: mean over vertices within searchlight
                    currentData = mean(currentData, 1); % (1, timePoints, conditions, sessions)
                    currentData = squeeze(currentData); % (timePionts, conditions, sessions)
                case 'spatiotemporal'
                    % Spatiotemporal patterns: all the data concatenated
                    currentData = reshape(currentData, [], size(currentData, 3), size(currentData, 4)); % (dataPoints, conditions, sessions)
            end%switch:userOptions.sensorSearchlightPatterns
            
            for session = 1:userOptions.nSessions
                
                searchlightRDM = searchlightRDM + squareform(pdist(squeeze(currentData(:,:,session))',userOptions.distance));
                
            end%for:sessions
            
        else % data regularization based on algorithm by Diedrichson et al 2011 - updated 12-12 IZ
            tempMesh = reshape(currentData, [], size(currentData, 3), size(currentData, 4));
            currentData = zeros(size(tempMesh, 1), size(tempMesh, 2) * size(tempMesh, 3)); % (data, conditions, sessions)
            
            % combining session-wise trials
            kk = 1;
            for j = 1:size(tempMesh,2)
                for i = 1:userOptions.nSessions
                    currentData(:, kk) = (tempMesh(:, j, i));
                    kk = kk + 1;
                end
            end
            
            r_matrix = g_matrix(zscore(squeeze(currentData(:,:)))', userOptions.nConditions, size(currentTimeWindow,2));
            searchlightRDM = searchlightRDM + (1 - r_matrix);
            
            if isnan(searchlightRDM) % sessions and conditions should be optimal
                error('Cannot calculate covariance matrix. Try reducing number of conditions');
            end
        end
        
        searchlightRDM = searchlightRDM / userOptions.nSessions;
        
        searchlightRDM = vectorizeRDM(searchlightRDM);
        
        % Locally store the full brain's worth of indexed RDMs. (just
        % lower triangle) added by IZ 09-12
        if strcmp(userOptions.groupStats, 'FFX')
%             searchlightRDMs(:,:,vertex, t) = single(tril(squareform(searchlightRDM)));
              searchlightRDMs.(['v_' num2str(vertex)]).(['t_' num2str(t)]).RDM = searchlightRDM;
        else
            searchlightRDMs = nan(1);
        end
        
        % TODO: Refactor this into general method so it can be used
        % TODO: anywhere
        if strcmpi(userOptions.distanceMeasure, 'Kendall_taua')
            rs = rankCorr_Kendall_taua(searchlightRDM', modelRDM_utv');
            % We are currently no using the p values created, though we are
            % saving them for some reason.  Anyway, there is currently no
            % way to get them for Kendall's tau_a, and until this becomes
            % requiresd, I'm not going to fix it :) - C
            ps = NaN;
        elseif userOptions.partial_correlation
            % TODO: Consider partialcorr with kendall's tau
            [rs, ps] = partialcorr(searchlightRDM', modelRDM_utv', control_for_modelRDMs', 'type', userOptions.distanceMeasure, 'rows','pairwise');
        else
            [rs, ps] = corr(searchlightRDM', modelRDM_utv', 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
        end
        
        smm_ps(vertex, t, :) = ps;
        smm_rs(vertex, t, :) = rs;
        
    end%for:t
    
    % Indicate progress every once in a while...
    if mod(vertex, floor(length(vertices) / 20)) == 0, fprintf('.'); end%if
    
end%for:vertices

if userOptions.fisher
    smm_rs = fisherTransform(smm_rs);
end%if

end%function
