% [glmMeshPaths, lagSTCMetadata] = ...
%     searchlightGLM(averageRDMPaths, models, dataSTCMetadata, userOptions ...
%                   ['lag', <lag_in_ms>])
%
% models: Is a nTimepoints x nModels struct with field .RDM
%
% dataSTCMetadata: Contains info about timing and vertices for the data, 
%                  it's necessary for applying appropriate lags to the
%                  models.
%
% lag: The lag offset for the model time courses in ms. Must be
%      non-negative.
%
% Based on scripts written by Li Su and Isma Zulfiqar.
%
% Cai Wingfield 2015-03 -- 2015-04
function [glmMeshPaths, lagSTCMetadatas] = searchlight_dynamicGLM_source(averageRDMPaths, models, slSTCMetadatas, userOptions, varargin)

    import rsa.*
    import rsa.rdm.*
    import rsa.util.*
    
    %% Parse inputs
    
    % 'lag'
    nameLag = 'lag';
    checkLag = @(x) (isnumeric(x) && (x >= 0));
    defaultLag = 0;
    
    % Set up parser
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.StructExpand  = false;
    
    % Parameters
    addParameter(ip, nameLag, defaultLag, checkLag);
    
    % Parse the inputs
    parse(ip, varargin{:});
    
    % Get some nicer variable names
    % The lag in ms
    lag_in_ms = ip.Results.lag; % 111
    
    [nTimepoints_models, nModels] = size(models);
    
    
    %% Begin
    
    for chi = 'LR'
    
    
        %% Prepare lag for the models

        prints('Computing appropriate lag for dynamic model GLM...');

        % The models are assumed to have the same number of timepoints as the
        % data, and the timepoints are assumed to be corresponding.

        % The timepoints in the model timelines and the timepoints in the data
        % timelines are assumed to be corresponding at 0 lag, though the models
        % will be  offset by the specified lag.

        % Remember that STCmetadata.tstep measures lag in SECONDS!
        timestep_in_ms = slSTCMetadatas.(chi).tstep * 1000;

        % Check if this lag is doable
        if mod(lag_in_ms, timestep_in_ms) ~= 0
            warns('The requested lag of %dms cannot be achieved, as the timestep is %dms.', lag_in_ms, timestep_in_ms);

            % If it's not achievable, we adjust it until it is
            desired_lag_in_steps = lag_in_ms / timestep_in_ms;
            % TODO: this takes the floor, but should really take the nearest?
            achievable_lag_in_steps = floor(desired_lag_in_steps);
            achievable_lag_in_ms = achievable_lag_in_steps * timestep_in_ms;
            warns('Using a lag of %dms instead.', achievable_lag_in_ms);
            lag_in_ms = achievable_lag_in_ms;
        end

        lag_in_timepoints = lag_in_ms / timestep_in_ms;
    
    
        %% Prepare lag STC metadata

        lagSTCMetadatas.(chi).tstep = slSTCMetadatas.(chi).tstep;
        lagSTCMetadatas.(chi).vertices = slSTCMetadatas.(chi).vertices;
        lagSTCMetadatas.(chi).tmax = slSTCMetadatas.(chi).tmax;
        lagSTCMetadatas.(chi).tmin = slSTCMetadatas.(chi).tmin + (lagSTCMetadatas.(chi).tstep * lag_in_timepoints);
        
        prints('Loading average RDM mesh from "%s"...', averageRDMPaths.(chi));
        
        average_slRDMs = directLoad(averageRDMPaths.(chi), 'average_slRDMs');
        
        prints('Applying lag to dynamic model timelines...');
    
        [nVertices, nTimepoints_data] = size(average_slRDMs);
        [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data);
    
        prints('Working at a lag of %dms, which corresponds to %d timepoints at this resolution.', lag_in_ms, lag_in_timepoints);
        
        % Preallocate
        glm_mesh(1:nVertices, 1:nTimepoints_overlap) = struct('betas', nan, 'deviance', nan, 'maxBeta', nan, 'maxBeta_i', nan);
        
        prints('Performing dynamic GLM in %sh hemisphere...', lower(chi));
        
        parfor t = 1:nTimepoints_overlap
            
            t_relative_to_data = t + lag_in_timepoints;
    
            % Temporarily dissable this warning
            warning_id = 'stats:glmfit:IllConditioned';
            warning('off', warning_id);

            prints('Working on timepoint %d/%d...', t, nTimepoints_overlap);
            
            for v = 1:nVertices
            
                % Fit the GLM at this point
                % TODO: In case the models are all zeros, this will merrily
                % TODO: produce meaningless betas along with a warning.
                % TODO: We should probably check for this first.
                [ ...
                      glm_mesh(v, t).betas ...
                    , glm_mesh(v, t).deviance ...
                    ...% TODO: do we need these stats? It sure is a huge 
                    ...% TODO: amount of data to keep in memory and read/
                    ...% TODO: write from disk
                    ...%, glm_mesh(v, t).stats ...
                    ] = glmfit( ...
                        modelStack{t}', ...
                        average_slRDMs(v, t_relative_to_data).RDM', ...
                        ...% TODO: Why are we making this assumption?
                        ...% TODO: What are the implications of this?
                        'normal');
                
                % TODO: In case of a tie, this takes the first beta.
                % TODO: It would be better to take a random one, perhaps
                % TODO: using rsa.util.chooseRandom() somehow.
                % TODO: Make it clear that the i-s are of the original
                % TODO: betas, not of the betas in the list (which also has
                % TODO: the first one as an all-ones beta).
                [glm_mesh(v, t).maxBeta, glm_mesh(v, t).maxBeta_i] = max(glm_mesh(v, t).betas(2:end));
                
            end%for:v
            
            % Re-enable warning
            warning('on', warning_id);
            
        end%for:t

        %% Save results
        
        glmMeshDir = fullfile(userOptions.rootPath, 'Meshes');
        glmMeshFilename = ['GLM_mesh_', lower(chi), 'h.mat'];
        glmMeshPaths.(chi) = fullfile(glmMeshDir, glmMeshFilename);
        
        prints('Saving GLM results for %sh hemisphere to %s...', lower(chi), glmMeshPaths.(chi));
        
        gotoDir(glmMeshDir);
        save('-v7.3', glmMeshPaths.(chi), 'glm_mesh');
        
        
        %% Save STCs
        
        % Preallocate
        max_beta_values = zeros(nVertices, nTimepoints_overlap);
        max_beta_model_is = zeros(nVertices, nTimepoints_overlap);
        deviances = zeros(nVertices, nTimepoints_overlap);
        all_data = zeros(nVertices, nTimepoints_overlap, numel(glm_mesh(1,1).betas));
        
        % Reshape data into interesting formats
        parfor t = 1:nTimepoints_overlap
            for v = 1:nVertices
                % We use 2:end because the GLM automatically puts an all-1s
                % vector in as a predictor, but we're not interested in
                % that.
                [max_beta_values(v, t), max_beta_model_is(v, t)] = max(glm_mesh(v, t).betas(2:end));
                deviances(v, t) = glm_mesh(v, t).deviance;
                all_data(v, t, :) = glm_mesh(v, t).betas;
            end
        end
        
        % Copy the lag-adjusted metadata structure.
        stc_metadata_to_save = lagSTCMetadatas.(chi);
        
        gotoDir(glmMeshDir);
        
        % Save max beta values at each vertex.
        stc_file_name = sprintf('max_betas-%sh.stc', lower(chi));
        prints('Saving maximum beta values to %s...', stc_file_name);
        stc_metadata_to_save.data = max_beta_values;
        mne_write_stc_file1(fullfile(glmMeshDir, stc_file_name), stc_metadata_to_save);
        
        % Save the index of the model with the max beta value at each
        % vertex (not including the all-1s predictor of glmfit).
        stc_file_name = sprintf('max_beta_model_is-%sh.stc', lower(chi));
        prints('Saving peak model indices to %s...', stc_file_name);
        stc_metadata_to_save.data = max_beta_model_is;
        mne_write_stc_file1(fullfile(glmMeshDir, stc_file_name), stc_metadata_to_save);
        
        % Save the devainces at each vertex.
        stc_file_name = sprintf('deviances-%sh.stc', lower(chi));
        prints('Saving deviances to %s...', stc_file_name);
        stc_metadata_to_save.data = deviances;
        mne_write_stc_file1(fullfile(glmMeshDir, stc_file_name), stc_metadata_to_save);
        
        % Save the beta values for each individual model at each vertex.
        prints('Saving individual model betas...');
        parfor model_i = 1:nModels
            stc_file_name = sprintf('model_%d_betas-%sh.stc', model_i, lower(chi));
            prints('Saving betas for model %d in %s...', model_i, stc_file_name);
            
            % We use model_i+1 here as the model indices are offset by 1 by
            % glmfit's all-1s vector.
            stc_metadata_to_save = all_data(:, :, model_i+1);
            mne_write_stc_file1(fullfile(glmMeshdir, stc_file_name), stc_metadata_to_save);
        end
        
    end%for:chi
    
end%function


% [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data)
%
% Given some dynamic models, and a specified offset in timepoints, produces
% a time-indexed cell array of model stacks suitable for glmfit.
%
% Cai Wingfield 2015-04
function [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data)

    import rsa.*
    import rsa.rdm.*

    [nTimepoints_models, nModels] = size(models);
    
    % We only look at timepoints where the data's timeline overlaps with
    % the models' lag-offset timelines.
    %
    %    (lag)>  |--------------------| lag-offset models
    % |--------------------| data
    %            .         .
    %            |---------|
    %                 ^
    %             only look
    %             in overlap
    nTimepoints_overlap = nTimepoints_data - lag_in_timepoints;
    
    model_size = size(models(1,1).RDM);
    
    % Make sure we're using ltv form.
    model_size = size(vectorizeRDM(zeros(model_size)));
    
    % Now at each timepoint we stack the models into a predictor matrix for
    % the GLM.
    % We are only looking in the first bit of the models' timelines, in the
    % places where there is also data.
    for t = 1:nTimepoints_overlap
        for model_i = 1:nModels
            modelStack{t}(model_i, :) = vectorizeRDM(models(t, model_i).RDM);
        end%for:model
    end%for:t
end%function
