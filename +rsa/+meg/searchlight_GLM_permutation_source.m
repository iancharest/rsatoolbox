% searchlight_thresholdGLM_source(averageRDMPaths, glm_paths, models, slSTCMetadatas, lagSTCMetadatas, nPermutations)
%
% Cai Wingfield 2015-04
function [p_paths, p_median_paths] = searchlight_GLM_permutation_source(averageRDMPaths, glm_paths, models, slSTCMetadatas, lagSTCMetadatas, nPermutations)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.stat.*
    import rsa.util.*
    
    
    %% Precompute permutations for speed
    
    % Indices for lower-triangular-form of RDM.
    lt_indices = 1:numel(vectorizeRDM(models(1).RDM));
    
    % Indices for squareform of RDM.
    sf_indices = squareform(lt_indices);
    
    % Preallocate
    lt_index_permutations = nan(numel(lt_indices), nPermutations);
    
    prints('Precomputing index permutations');
    
    parfor p = 1:nPermutations
        lt_index_permutations(:, p) = squareform(randomizeSimMat(sf_indices));
    end
    
    
    %% Both hemispheres separately.
    
    for chi = 'LR'
        
        %% Load data RDMs
        
        prints('Loading %sh data RDMs from "%s"...', lower(chi), averageRDMPaths.(chi));
        average_slRDMs = directLoad(averageRDMPaths.(chi), 'average_slRDMs');
        
        [nVertices, nTimepoints_data] = size(average_slRDMs);
    
        lag_in_timepoints = (lagSTCMetadatas.(chi).tmin - slSTCMetadatas.(chi).tmin) / lagSTCMetadatas.(chi).tstep;

        [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data);

        nModels = numel(modelStack);
        % + 1 for that all-1s predictor
        nBetas = nModels + 1;

        
        %% Calculate pooled-over-time distributions of betas at each vertex

        % We'll pool across timepoints. So we'll make h0_betas a
        % nModels-x-(nPermutations*nTimepoints_overlap)-sized matrix.
        % This is based on the assumption that the distributions of
        % beta values should be independent of time.
        % We may (or may not) want to make the same assumption about
        % space, but we won't do that for now.
        
        % Temporarily disable this warning
        warning_id = 'stats:glmfit:IllConditioned';
        warning('off', warning_id);

        % Preallocate
        h0_betas = zeros(nVertices, nTimepoints_overlap, nBetas, nPermutations);

        prints('Computing beta null distrubitions at %d vertices, %d permutations each...', nVertices, nPermutations);
            
        parfor t = 1:nTimepoints_overlap
            t_relative_to_data = t + lag_in_timepoints;

            for p = 1:nPermutations
                for v = 1:nVertices

                    scrambled_data_rdm = average_slRDMs(v, t_relative_to_data).RDM(lt_index_permutations(:, p));

                    h0_betas(v, t, :, p) = glmfit( ...
                        modelStack{t}', ...
                        scrambled_data_rdm', ...
                        'normal');
                end%for
            end%for
        end%for
        
        % Re-enable warning
        warning('on', warning_id);
        
        % Reshape h0-betas
        % This will now be nVertices x pooled-values sized
        h0_betas = reshape(h0_betas, nVertices, nTimepoints_overlap * nPermutations * nBetas);


        %% Calculate p values
        
        prints('Calculating p-values at each vertex...');
        
        %% Load existing data
        
        prints('Loading actual %sh beta values from "%s"...', lower(chi), glmMeshPaths.(chi));
        
        glm_mesh_betas = directLoad([glm_paths.betas.(chi) '.mat'], 'glm_mesh_betas');
        glm_mesh_betas_median = directLoad([glm_paths.betas_median.(chi) '.mat'], 'glm_mesh_betas_median');

        % Preallocate
        p_mesh = ones(nVetices, nTimepoints_overlap, nModels);
        p_mesh_median = ones(nVertices, nModels);

        parfor m = 1:nModels
            prints('Computing p-values for model %d of %d...', m, nModels);
            p_mesh_slice = ones(nVertices, nTimepoints_overlap);
            p_mesh_median_slice = ones(nVertices, 1);
            for v = 1:nVertices
                for t = 1:nTimepoints_overlap
                    % +1 because of that annoying forced all-1s model.
                    p_mesh_slice(v, t) = 1 - portion(h0_betas(m, :), glm_mesh_betas(v, t, m + 1));
                end
                p_mesh_median_slice(v) = 1 - portion(h0_betas(m, :), glm_mesh_betas_median(v, m + 1));
            end
            p_mesh(:, :, m) = p_mesh_slice; %#ok<PFOUS> it's saved
            p_mesh_median(:, m) = p_mesh_median_slice; %#ok<PFOUS> it's saved
        end

        % Save results
        p_file_name = sprintf('p_mesh-%sh', lower(chi));
        p_paths.(chi) = fullfile(glmMeshDir, p_file_name);
        
        p_median_file_name = sprintf('p_mesh_median-%sh', lower(chi));
        p_median_paths.(chi) = fullfile(glmMeshDir, p_median_file_name);

        prints('Saving p-meshes to "%s"...', glmMeshDir);
        save(p_paths.(chi), 'p_mesh');
        save(p_median_paths.(chi), 'p_mesh_median');
    
    end%for:chi

end%function
