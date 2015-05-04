% [thresholded_glm_paths] = searchlight_GLM_threshold_source(glm_paths, h0_pooled_paths, lagSTCMetadatas, threshold)
%
% Cai Wingfield 2015-04 -- 2015-05
function [thresholded_glm_paths] = searchlight_GLM_threshold_source(glm_paths, h0_pooled_paths, lagSTCMetadatas, threshold, userOptions)
    
    import rsa.*
    import rsa.util.*
    import rsa.meg.*

    for chi = 'LR'
    
        thresholded_glm_paths.max_beta_is.(chi) = [glm_paths.max_beta_is.(chi) '_thresholded-' lower(chi) 'h'];
        thresholded_glm_paths.betas_model.(chi) = [glm_paths.betas_model.(chi) '_thresholded-' lower(chi) 'h'];
        
        %% Load null distributions
        
        prints('Loading %sh null distributions...', lower(chi));
        
        h0_betas = directLoad(h0_pooled_paths.(chi));
        
        [nVertices, poolSize] = size(h0_betas);
        
        %% Get threshold betas
        
        prints('Selecting top %2.1f centile of the null distribution at each vertex.', 100 * (1 - threshold));

        beta_thresholds = zeros(nVertices, 1);
        for v = 1:nVertices
           beta_thresholds(v) = quantile(h0_betas(v, :), 1 - threshold);
        end
        
        prints('This gives a median GLM coefficient threshold of %f over vertices.', median(beta_thresholds));

        
        %% Load beta values
        
        prints('Loading actual %sh beta values...', lower(chi));

        % (vertices, timepoints, models+1)
        glm_mesh_betas       = directLoad(glm_paths.betas.(chi),       'glm_mesh_betas');
        % (vertices, timepoints)
        glm_mesh_max_beta_is = directLoad(glm_paths.max_beta_is.(chi), 'glm_mesh_max_beta_is');
        
        [nVertices, nTimepoints_overlap, nBetas] = size(glm_mesh_betas);
        nModels = nBetas - 1;
        
        
        %% Threshold betas
        
        prints('Thresholding %sh beta meshes...', lower(chi));
        
        for v = 1:nVertices
            values_this_vertex = glm_mesh_betas(v, :, :);
            
            % Threshold max beta is
            for t = 1:nTimepoints_overlap
                max_beta_i = glm_mesh_max_beta_is(v, t);
                max_beta_value = values_this_vertex(1, t, max_beta_i + 1);
                if max_beta_value <= beta_thresholds(v)
                    max_beta_i = 0;
                end
                glm_mesh_max_beta_is(v, t) = max_beta_i;
            end
            
            % Threshold beta values
            values_this_vertex(values_this_vertex <= beta_thresholds(v)) = 0;
            glm_mesh_betas(v, :, :) = values_this_vertex;
        end
       
        
        %% Save thresholded meshes
        
        prints('Saving thresholded %sh values...', lower(chi));
        
        % Individual models
        for m = 1:nModels
            write_stc_file( ...
                lagSTCMetadatas.(chi), ...
                glm_mesh_betas(:, :, m + 1), ...
                sprintf([glm_paths.betas_model.(chi) '_thresholded-' lower(chi) 'h.stc'], m));
        end
        % Max betas
        write_stc_file( ...
            lagSTCMetadatas.(chi), ...
            glm_mesh_max_beta_is, ...
            [thresholded_glm_paths.max_beta_is.(chi) '.stc']);
        
    end%for:chi
end%function
