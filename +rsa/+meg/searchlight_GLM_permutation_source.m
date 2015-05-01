% [p_paths, p_median_paths] = searchlight_GLM_permutation_source(RDMPaths, glm_paths, models, slSTCMetadatas, lagSTCMetadatas, nPermutations, userOptions)
%
% Cai Wingfield 2015-04
function [p_paths, p_median_paths] = searchlight_GLM_permutation_source(RDMPaths, glm_paths, models, slSTCMetadatas, lagSTCMetadatas, nPermutations, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.stat.*
    import rsa.util.*
    
    
    %% Things to be returned whether or not work is done
    
    for chi = 'LR'
        
        glmMeshDir = fullfile(userOptions.rootPath, 'Meshes');
        
        p_file_name = sprintf('p_mesh-%sh', lower(chi));
        p_paths.(chi) = fullfile(glmMeshDir, p_file_name);
        
        p_median_file_name = sprintf('p_mesh_median-%sh', lower(chi));
        p_median_paths.(chi) = fullfile(glmMeshDir, p_median_file_name);
        
    end%for
    
    
    %% Precompute permutations for speed
    
    prints('Precomputing RDM index permutations...');
    
    % Indices for lower-triangular-form of RDM.
    lt_indices = 1:numel(vectorizeRDM(models(1).RDM));
    
    % Indices for squareform of RDM.
    sf_indices = squareform(lt_indices);
    
    % Preallocate
    lt_index_permutations = nan(numel(lt_indices), nPermutations);
    
    % Generate some permutations
    for p = 1:nPermutations
        lt_index_permutations(:, p) = squareform(randomizeSimMat(sf_indices));
    end
    
    
    %% Both hemispheres separately.
    
    for chi = 'LR'
        
        %% Load data RDMs
        
        prints('Loading %sh data RDMs from "%s"...', lower(chi), RDMPaths.(chi));
        
        slRDMs = directLoad(RDMPaths.(chi));
        
        [nVertices, nTimepoints_data] = size(slRDMs);
        lag_in_timepoints = (lagSTCMetadatas.(chi).tmin - slSTCMetadatas.(chi).tmin) / lagSTCMetadatas.(chi).tstep;

        [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data);

        nModels = size(modelStack{1}, 1);
        % + 1 for that all-1s predictor
        nBetas = nModels + 1;

        
        %% Calculate pooled-over-time distributions of betas at each vertex

        % Preallocate
        h0_betas = zeros(nVertices, nTimepoints_overlap, nBetas, nPermutations);

        prints('Computing beta null distributions at %d vertices...', nVertices);
            
        parfor t = 1:nTimepoints_overlap
        
            % Temporarily disable this warning
            %warning_id = 'stats:glmfit:IllConditioned';
            w = warning('off', 'all');
            
            prints('Timepoint %d of %d...', t, nTimepoints_overlap);
            
            t_relative_to_data = t + lag_in_timepoints;
            
            for v = 1:nVertices
                
                unscrambled_data_rdm = slRDMs(v, t_relative_to_data).RDM;

                for p = 1:nPermutations
		
                    scrambled_data_rdm = unscrambled_data_rdm(lt_index_permutations(:, p));

                    h0_betas(v, t, :, p) = glmfit( ...
                        modelStack{t}', ...
                        scrambled_data_rdm', ...
                        'normal');
                end%for
                
                % Occasional feedback
                if feedback_throttle(10, v, nVertices)
                    prints('%2.0f%% of vertices covered for timepoint %d.', percent(v, nVertices), t);
                end
            end%for
        
            % Re-enable warning
            %warning('on', warning_id);
            warning(w);
        end%for
        
        % Save null-distributions pre pooling
        gotoDir(userOptions.rootPath, 'Stats');
        save(sprintf('unpooled-h0-%sh', lower(chi)), 'h0_betas', '-v7.3');
        
        % We'll pool the distrubution across timepoints and permutations.
        % This is based on the assumption that the distributions of
        % beta values should be independent of time.
        % We may (or may not) want to make the same assumption about
        % space, but we won't do that for now.
        
        % But first we need to slice out the betas for the all-1s
        % predictor.
        h0_betas = h0_betas(:, :, 2:end, :);
        
        % We want the distribution of maximum-over-models betas at each
        % vertex.
        % (nVertices, nTimepoints_overlap, nPermutations)
        h0_betas = max(h0_betas, 3);
        % (nVertices, nTimepoints_overlap * nPermutations)
        h0_betas = reshape(h0_betas, ...
            nVertices, nTimepoints_overlap * nPermutations);
        
        % Save null-distributions post pooling
        save(sprintf('pooled-h0-%sh', lower(chi)), 'h0_betas', '-v7.3');

        
        %% Load existing data
        
        prints('Loading actual %sh beta values...', lower(chi));
        
        glm_mesh_betas = directLoad([glm_paths.betas.(chi) '.mat'], 'glm_mesh_betas');
        glm_mesh_betas_median = directLoad([glm_paths.betas_median.(chi) '.mat'], 'glm_mesh_betas_median');
        
        
        %% Calculate p values
        
        prints('Calculating p-values at each vertex...');

        % Preallocate
        p_mesh = ones(nVertices, nTimepoints_overlap, nModels);
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
        
        %% Save results

        prints('Saving p-meshes...');
        save(p_paths.(chi), 'p_mesh');
        save(p_median_paths.(chi), 'p_mesh_median');
    
    end%for:chi

end%function
