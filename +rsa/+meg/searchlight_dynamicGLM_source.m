% [glm_paths, lagSTCMetadata] = ...
%     searchlightGLM(RDMPaths, models, dataSTCMetadata, userOptions ...
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
% Cai Wingfield 2015-03 -- 2015-06
function [glm_paths, lagSTCMetadatas] = searchlight_dynamicGLM_source(RDMPaths, models, slSTCMetadatas, userOptions, varargin)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.stat.*
    import rsa.util.*
    
    %% Parse inputs
    
    % 'lag'
    nameLag = 'lag';
    checkLag = @(x) (isnumeric(x) && (x >= 0));
    defaultLag = 0;
    
    % 'file-prefix'
    nameFilePrefix = 'fileprefix';
    checkFilePrefix = @(x) (ischar(x));
    defaultFilePrefix = '';
    
    % Set up parser
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.StructExpand  = false;
    
    % Parameters
    addParameter(ip, nameLag, defaultLag, checkLag);
    addParameter(ip, nameFilePrefix, defaultFilePrefix, checkFilePrefix);
    
    % Parse the inputs
    parse(ip, varargin{:});
    
    % Get some nicer variable names
    
    % The lag in ms
    lag_in_ms = ip.Results.(nameLag);
    
    % The file name prefix
    file_name_prefix = ip.Results.(nameFilePrefix);
    
    
    %% Set up values to be returned, whether or not any work is really done.
    
    for chi = 'LR'
        
        %% Where to save results
        
        % Directory
        glmMeshDir = fullfile(userOptions.rootPath, 'Meshes');
        gotoDir(glmMeshDir);
        
        % The same paths will be used for mat files and stc files, the only
        % differences being the extension.
        
        % Paths
        glm_paths.betas.(chi) = fullfile(glmMeshDir, ...
            [file_name_prefix 'GLM_mesh_betas-',       lower(chi), 'h']);
        glm_paths.deviances.(chi) = fullfile(glmMeshDir, ...
            [file_name_prefix 'GLM_mesh_deviances-',   lower(chi), 'h']);
        
        % Paths model template
        glm_paths.betas_model.(chi) = fullfile(glmMeshDir, ...
            [file_name_prefix 'GLM_mesh_betas_model_%d-', lower(chi), 'h']);
        
        
        %% Prepare lag for the models

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
        % tmin is increased by...
        lagSTCMetadatas.(chi).tmin = slSTCMetadatas.(chi).tmin + ...
            ...% timesteps equal to...
            (lagSTCMetadatas.(chi).tstep * ( ...
                ...% the fixed lag we apply.
                lag_in_timepoints));
        
    end
    
    
    %% Check for overwrites
    
    promptOptions.functionCaller = 'searchlight_dynamicGLM_source';
    file_i = 1;
    for chi = 'LR'
        fileNames = fieldnames(glm_paths);
        for file_name_i = 1:numel(fileNames)
            fileName = fileNames{file_name_i};
            promptOptions.checkFiles(file_i).address = [glm_paths.(fileName).(chi) '.mat'];
            file_i = file_i + 1;
        end
    end
    % Some of the file names are templates, so we don't want these 'files'
    % failing to be detected to result in rerunning in every case.
    promptOptions.quantification = 'existential';
    promptOptions.defaultResponse = 'R';
    
    overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
    
    %% Begin

    if overwriteFlag
    
        [nTimepoints_models, nModels] = size(models);

        for chi = 'LR'

            prints('Loading RDM mesh from "%s"...', RDMPaths.(chi));

            slRDMs = directLoad(RDMPaths.(chi));

            prints('Applying lag to dynamic model timelines...');

            [nVertices, nTimepoints_data] = size(slRDMs);
            [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data);

            prints('Working at a lag of %dms, which corresponds to %d timepoints at this resolution.', lag_in_ms, lag_in_timepoints);

            % Preallocate.
            glm_mesh_betas = nan(nVertices, nTimepoints_overlap, nModels + 1);
            glm_mesh_deviances = nan(nVertices, nTimepoints_overlap);

            % Tell the user what's going on.
            prints('Performing dynamic GLM in %sh hemisphere...', lower(chi));

            parfor t = 1:nTimepoints_overlap

                prints('Working on timepoint %d/%d...', t, nTimepoints_overlap);

                % Temporarily disable this warning
                warning_id = 'stats:glmfit:IllConditioned';
                warning('off', warning_id);

                % The timelines for the data and the models are offset
                t_relative_to_data = t ...
                    ...% by a fixed lag
                    + lag_in_timepoints;

                for v = 1:nVertices
                    [glm_mesh_betas(v, t, :), glm_mesh_deviances(v, t)] = ...
                        glmfit( ...
                            modelStack{t}', ...
                            slRDMs(v, t_relative_to_data).RDM', ...
                            'normal');
                end%for:v

                % Re-enable warning
                warning('on', warning_id);

            end%for:t


            %% Save results
            prints('Saving GLM results for %sh hemisphere to "%s"...', lower(chi), glmMeshDir);

            %% Save full results
            save('-v7.3', glm_paths.betas.(chi), 'glm_mesh_betas');

            %% Save per-model results
            for model_i = 1:nModels
                write_stc_file( ...
                    lagSTCMetadatas.(chi), ...
                    squeeze(glm_mesh_betas(:, :, model_i + 1)), ...
                    [sprintf(glm_paths.betas_model.(chi), model_i) '.stc']);
            end

            %% Save summary results
            save('-v7.3', glm_paths.deviances.(chi), 'glm_mesh_deviances');
            write_stc_file( ...
                lagSTCMetadatas.(chi), ...
                glm_mesh_deviances, ...
                [glm_paths.deviances.(chi) '.stc']);

        end%for:chi
    end%if
end%function
