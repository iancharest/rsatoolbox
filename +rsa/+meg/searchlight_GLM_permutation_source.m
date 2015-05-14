% [h0_paths] = searchlight_GLM_permutation_source(RDMPaths, models, slSTCMetadatas, lagSTCMetadatas, nPermutations, userOptions)
%
% Cai Wingfield 2015-04
function [h0_paths] = searchlight_GLM_permutation_source(RDMPaths, models, slSTCMetadatas, lagSTCMetadatas, nPermutations, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.stat.*
    import rsa.util.*
    
    %% Things to be returned whether or not work is done
    
    for chi = 'LR'
        % Where to save results
        glmStatsDir = fullfile(userOptions.rootPath, 'Stats');
        h0_paths.(chi) = fullfile(glmStatsDir, sprintf('h0-%sh.mat', lower(chi)));
    end%for
    
    
    %% Check for overwrites
    
    promptOptions.functionCaller = 'searchlight_GLM_permutation_source';
    file_i = 1;
    for chi = 'LR'
        promptOptions.checkFiles(file_i).address = h0_paths.(chi);
        file_i = file_i + 1;
    end
    promptOptions.defaultResponse = 'R';
    
    overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
    
    %% Begin
    
    if overwriteFlag

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


            %% Calculate distributions of betas at each vertex

            % Preallocate
            h0_betas = zeros( ...
                nVertices, ...
                nTimepoints_overlap, ...
                ...% -1 for the all-1s beta which we will ignore
                nBetas - 1, ...
                nPermutations);

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
                        
                        betas = glmfit( ...
                            modelStack{t}', ...
                            scrambled_data_rdm', ...
                            'normal');

                        h0_betas(v, t, :, p) = betas(2:end); %#ok<PFOUS> it's saved
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


            %% Pool and save H0-distributions

            % Ssave null-distributions pre pooling
            gotoDir(userOptions.rootPath, 'Stats');
            save(h0_paths.(chi), 'h0_betas', '-v7.3');

        end%for:chi
    end
end%function
