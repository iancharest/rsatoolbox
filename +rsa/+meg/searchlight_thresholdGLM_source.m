function searchlight_thresholdGLM_source(averageRDMPaths, glmMeshPaths, models, slSTCMetadatas, lagSTCMetadatas, userOptions)

    import rsa.*
    import rsa.rdm.*
    import rsa.util.*

    % Use L only
    lag_in_timepoints = (lagSTCMetadatas.L.tmin - slSTCMetadatas.L.tmin) / lagSTCMetadatas.L.tstep;
    nTimepoints_data = (slSTCMetadatas.L.tmax - slSTCMetadatas.L.tmin) - slSTCMetadatas.L.tstep;

    [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data);
    
    nModels = numel(modelStack);

    %% Calculate pooled-over-time distributions of betas at each vertex

        % We'll pool across timepoints. So we'll make h0_betas a
        % nModels-x-(nPermutations*nTimepoints_overlap)-sized matrix.
        % This is based on the assumption that the distributions of
        % beta values should be independnt of time.
        % We may or may not want to make the same assumption about
        % space, but we won't do that for now.

        prints('Computing beta distrubitions over %d vertices...', nVertices);
        parfor v = 1:nVertices

            % Calculate beta-distributions

            % We're storing timepoints and permutations and 
            h0_i = 1;
            for p = 1:nPermutations
                prints('Permutation %d of %d...', p, nPermutations);

                for t = 1:nTimepoints_overlap
                    t_relative_to_data = t + lag_in_timepoints;

                    permuted_data_rdm = randomizeSimMat(average_slRDMs(v, t_relative_to_data).RDM);

                    h0_betas(:, h0_i) = glmfit( ... % (nModels, nPermutations)
                        modelStack{t}', ...
                        permuted_data_rdm', ...
                    ...% TODO: Why are we making this assumption?
                    ...% TODO: What are the implications of this?
                        'normal');
                    h0_i = h0_i + 1;
                end
            end
        end

        % Calculate p-vales
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
        save(p_path, 'p_mesh');

    end

end%function
