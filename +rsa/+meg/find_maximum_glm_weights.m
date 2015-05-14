% Cai Wingfield 2015-05
function [glm_paths] = find_maximum_glm_weights(glm_paths, lagSTCMetadatas, userOptions)
    
    import rsa.*
    import rsa.util.*
    import rsa.meg.*
    
    glmMeshDir = fullfile(userOptions.rootPath, 'Meshes');

    for chi = 'LR'
        
        glm_paths.max_betas.(chi) = fullfile(glmMeshDir, ...
            ['GLM_mesh_max_betas-',   lower(chi), 'h']);
        glm_paths.max_beta_is.(chi) = fullfile(glmMeshDir, ...
            ['GLM_mesh_max_beta_is-', lower(chi), 'h']);
        
        %% Load beta values
        
        prints('Loading %sh beta values...', lower(chi));

        % (vertices, timepoints, models+1)
        glm_mesh_betas = directLoad(glm_paths.betas.(chi), 'glm_mesh_betas');
        
        [nVertices, nTimepoints_overlap, nBetas] = size(glm_mesh_betas);
        
        %% Computing max betas
        
        [glm_mesh_max_betas, glm_mesh_max_beta_is] = max(glm_mesh_betas(:, :, 2:end), [], 3);

        
        %% Saving max betas
        
        % Matlab format
        save('-v7.3', glm_paths.max_betas.(chi),   'glm_mesh_max_betas');
        save('-v7.3', glm_paths.max_beta_is.(chi), 'glm_mesh_max_beta_is');

        % STC format
        write_stc_file( ...
            lagSTCMetadatas.(chi), ...
            glm_mesh_max_betas, ...
            [glm_paths.max_betas.(chi) '.stc']);
        write_stc_file( ...
            lagSTCMetadatas.(chi), ...
            glm_mesh_max_beta_is, ...
            [glm_paths.max_beta_is.(chi) '.stc']);

        
    end%for:chi
end%function
