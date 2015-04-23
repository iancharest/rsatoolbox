function searchlight_thresholdGLM_source(averageRDMPaths, glmMeshPaths, models, slSTCMetadatas, lagSTCMetadatas, userOptions)

    import rsa.*
    import rsa.rdm.*
    import rsa.util.*
    
    % Use L only
    lag_in_timepoints = (lagSTCMetadatas.L.tmin - slSTCMetadatas.L.tmin) / lagSTCMetadatas.L.tstep;
    nTimepoints_data = (slSTCMetadatas.L.tmax - slSTCMetadatas.L.tmin) - slSTCMetadatas.L.tstep;
    
    [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data);
    
    nModels = numel(modelStack);
    
    for chi = 'LR'
        
        %% Load data RDMs
        
        prints('Loading %sh data RDMs from "%s"...', lower(chi), averageRDMPaths.(chi));
        average_slRDMs = directLoad(averageRDMPaths.(chi), 'average_slRDMs');

        %% Calculate pooled-over-time distributions of betas at each vertex

        % We'll pool across timepoints. So we'll make h0_betas a
        % nModels-x-(nPermutations*nTimepoints_overlap)-sized matrix.
        % This is based on the assumption that the distributions of
        % beta values should be independnt of time.
        % We may (or may not) want to make the same assumption about
        % space, but we won't do that for now.
        
        % Temporarily disable this warning
        warning_id = 'stats:glmfit:IllConditioned';
        warning('off', warning_id);

        pooled_dist_size = nPermutations * nTimepoints_overlap * nModels;

        prints('Computing beta null distrubitions at %d vertices, %d permutations each...', nVertices, nPermutations);
        prints('(This will give the distribution at each vertex a total of %d values.)', pooled_dist_size);

        % Preallocate
        h0_betas = zeros(nVertices, pooled_dist_size);

        for v = 1:nVertices

            h0_betas_this_vertex = zeros(nPermutations, nTimepoints_overlap * nModels);

            parfor p = 1:nPermutations
                prints('Permutation %d of %d...', p, nPermutations);

                h0_betas_this_perm = zeros(nTimepoints_overlap, nModels);

                for t = 1:nTimepoints_overlap

                    t_relative_to_data = t + lag_in_timepoints;

                    scrambled_data_rdm = randomizeSimMat(average_slRDMs(v, t_relative_to_data).RDM);

                    h0_betas_this_perm(t, :) = glmfit( ...
                        modelStack{t}', ...
                        scrambled_data_rdm', ...
                        'normal');
                end%for:t

                h0_betas_this_vertex(p, :) = h0_betas_this_perm(:);

            end%for:p

            h0_betas(v, :) = h0_betas_this_vertex(:);

            % Give feedback every once in a while
            if mod(v,20) == 0
                prints('Null distibutions simulated for %d%% of vertices.', 100 * v / nVertices);
            end

        end%for:v
        
        % Re-enable warning
        warning('on', warning_id);


        %% Calculate p values
        
        prints('Calculating p-values at each vertex...');
        
        %% Load existing data
        
        prints('Loading actual %sh beta values from "%s"...', lower(chi), glmMeshPaths.(chi));
        glm_mesh = directLoad(glmMeshPaths.(chi), 'glm_mesh');

        % Preallocate
        p_mesh = ones(nVetices, nTimepoints_overlap, nModels);

        parfor m = 1:nModels
            prints('Computing p-values for model %d of %d...', m, nModels);
            for v = 1:nVertices
                for t = 1:nTimepoints_overlap
                    % +1 because of that annoying forced all-1s model.
                    p_mesh(v, t, m) = 1 - portion(h0_betas(m, :), glm_mesh(v, t).betas(m + 1));
                end
            end
        end

        % Save results
        p_file_name = sprintf('p_mesh-%sh', lower(chi));
        p_path = fullfile(glmMeshDir, p_file_name);

        prints('Saving p-meshes to "%s"...', p_path);
        save(p_path, 'p_mesh');
    
    end%for:chi

end%function
