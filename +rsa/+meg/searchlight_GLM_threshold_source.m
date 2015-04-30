function [thresholded_glm_paths] = searchlight_GLM_threshold_source(glm_paths, p_paths, p_median_paths, lagSTCMetadatas, threshold)
    
    import rsa.*
    import rsa.util.*
    import rsa.meg.*

    for chi = 'LR'
        
        % (vertices, timepoints, models)
        glm_mesh_max_beta_is = directLoad(glm_paths.max_beta_is.(chi), 'glm_mesh_max_beta_is');
        p_mesh = directLoad(p_paths.(chi), 'p_mesh');
        p_mesh = min(p_mesh, [], 3);
        
        % (vertices, models)
        glm_mesh_max_beta_is_median = directLoad(glm_paths.max_beta_is_median.(chi), 'glm_mesh_max_beta_is_median');
        p_mesh_median = directLoad(p_median_paths.(chi), 'p_mesh_median');
        p_mesh_median = min(p_mesh_median, [], 2);
        
        % p-threshold beta meshes.
        glm_mesh_max_beta_is(p_mesh >= threshold) = 0;
        glm_mesh_max_beta_is_median(p_mesh_median >= threshold) = 0;
        
        % Save them again
        write_stc_file( ...
            lagSTCMetadatas.(chi), ...
            glm_mesh_max_beta_is, ...
            [glm_paths.max_beta_is.(chi) '_thresholded-' lower(chi) 'h.stc']);
        write_stc_snapshot( ...
            lagSTCMetadatas.(chi), ...
            glm_mesh_max_beta_is_median, ...
            [glm_paths.max_beta_is_median.(chi) '_thresholded-' lower(chi) 'h.stc']);
        
    end%for:chi
end%function
