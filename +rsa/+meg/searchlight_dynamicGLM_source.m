% [glmMeshPaths] = searchlightGLM(averageRDMPaths, models, dataSTCMetadata, userOptions ...
%                                ['lag', <lag_in_ms>])
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
% Cai Wingfield 2015-03
function [glmMeshPaths] = searchlight_dynamicGLM_source(averageRDMPaths, models, dataSTCMetadata, userOptions, varargin)

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
    
    
    %% Prepare lag for the models
    
    prints('Computing appropriate lag for dynamic model GLM...');
    
    % The models are assumed to have the same number of timepoints as the
    % data, and the timepoints are assumed to be corresponding.  All
    % computations of lag are based on these assumptions, and won't work at
    % all if they are violated.
    
    timestep_in_ms = dataSTCMetadata.tstep;
    
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
    
    lag_in_timePonts = lag_in_ms / timestep_in_ms;
    
    
    %% Vectorise and apply lag to models
    
    prints('Applying lag to dynamic models...');
    
    [nTimepoints, nModels] = size(models);
    model_size = size(models(1,1).RDM);
    
    % Make sure we're using ltv form.
    model_size = size(vectorizeRDM(zeros(model_size)));
    
    % Preallocate the models as NaNs.
    % We'll then overwrite them as necessary
    offset_models(1:nTimepoints, 1:nModels) = struct('RDM', nan(model_size));
    
    for model_i = 1:nModels
        % This index tracks the timepoint indices of the data.
        % The model's timepoint indices will be offset from this.
       for datalocked_timepoint = 1:nTimepoints
           offset_timepoint = datalocked_timepoint + lag_in_timePonts;
           
           % We don't want to exceed the bounds, so we make the check each
           % time.
           if offset_timepoint <= nTimepoints
               offset_models(offset_timepoint, model_i).RDM = models(datalocked_timepoint, model_i).RDM;
           end
       end
    end
    
    % `offset_models` now has the the same models as `models`, but with the
    % time indices positively shifted by lag.  This means we may lose some
    % of the later models, but this is ok, as with the hypothesised lag,
    % they are not predicted to correlate with any data in our epoch.
    
    % Now at each timepoint we stack the models into a predictor matrix for
    % the GLM.
    for t = 1:nTimepoints
        for model_i = 1:nModels
            modelStack{t}(model_i, :) = vectorizeRDM(models(model_i).RDM);
        end%for:model
    end
    
    
    %% Begin
    
    for chi = 'LR'
        
        prints('Performing dynamic GLM in %sh hemisphere...', lower(chi));
        
        average_slRDMs = directLoad(averageRDMPaths.(chi), 'average_slRDMs');
        
        [nVertices, nTimepoints] = size(average_slRDMs);
        
        % preallocate
        glm_mesh(1:nVertices, 1:nTimepoints) = struct();
        
        parfor t = 1:nTimepoints
            prints('Working on timepoint %d/%d...', t, nTimepoints);
            
            for v = 1:nVertices
            
                % Fit the GLM at this point
                % TODO: In case the models are all zeros, this will merrily
                % TODO: produce meaningless betas along with a warning.
                % TODO: We should probably check for this first.
                [ ...
                    glm_mesh(v, t).betas, ...
                    glm_mesh(v, t).deviance, ...
                    glm_mesh(v, t).stats] = glmfit( ...
                        modelStack{t}', ...
                        average_slRDMs(v, t).RDM', ...
                        ...% TODO: Why are we making this assumption?
                        ...% TODO: What are the implications of this?
                        'normal');
                
                % TODO: In case of a tie, this takes the first beta.
                % TODO: It would be better to take a random one, perhaps
                % TODO: using rsa.util.chooseRandom() somehow.
                [glm_mesh(v, t).maxBeta, glm_mesh(v, t).maxBeta_i] = max(glm_mesh(v, t).betas);
                
            end%for:v
        end%for:t

        %% Save results
        
        glmMeshDir = fullfile(userOptions.rootPath, 'Meshes');
        glmMeshFilename = ['GLM_mesh_', lower(chi), 'h.mat'];
        glmMeshPaths.(chi) = fullfile(glmMeshDir, glmMeshFilename);
        
        prints('Saving GLM results for %sh hemisphere to %s...', lower(chi), glmMeshPaths.(chi));
        
        gotoDir(glmMeshDir);
        save('-v7.3', glmMeshPaths.(chi), 'glm_mesh');
        
    end%for:chi
end%function
