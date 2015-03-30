function [glmMeshPaths] = searchlightGLM(averageRDMPaths, models, userOptions)

    import rsa.*
    import rsa.rdm.*
    import rsa.util.*
    
    models = vectorizeRDMs(models);
    nModels = numel(models);
    
    modelStack = nan(nModels, numel(models(1).RDM));
    for model_i = 1:nModels
        modelStack(model_i, :) = models(model_i).RDM;
    end%for:model
    
    for chi = 'LR'
        average_slRDMs = directLoad(averageRDMPaths.(chi), 'average_slRDMs');
        
        [nVertices, nTimepoints] = size(average_slRDMs);
        
        % preallocate
        glm_mesh(1:nVertices, 1:nTimepoints) = struct();
        
        parfor t = 1:nTimepoints
            for v = 1:nVertices
            
                % Fit the GLM at this point
                [ ...
                    glm_mesh(v, t).betas, ...
                    glm_mesh(v, t).deviance, ...
                    glm_mesh(v, t).stats] = glmfit( ...
                        modelStack, ...
                        average_slRDMs(v, t).RDM, ...
                        ...% TODO: Why are we making this assumption?
                        'normal');
                
                % TODO: In case of a tie, this takes the first beta.
                % TODO: It would be better to take a random one, perhaps
                % TODO: using rsa.util.chooseRandom() somehow.
                [glm_mesh(v, t).maxBeta, glm_mesh(v, t).maxBeta_i] = max(glm_mesh(v, t).betas);
                
            end%for:v
        end%for:t

        glmMeshDir = fullfile(userOptions.rootPath, 'Meshes');
        glmMeshFilename = ['GLM_mesh_', lower(chi), 'h.mat'];
        glmMeshPaths.(chi) = fullfile(glmMeshDir, glmMeshFilename);
        gotoDir(glmMeshDir);
        save('-v7.3', glmMeshPaths.(chi), 'glm_mesh');
        
    end%for:chi
    
end%function
