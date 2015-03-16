% sliding time window analysis for ROIs
% Written by Isma Zulfiqar 12-12 -- Updated 03/13
% Updated by Fawad 03-2014, 10-2014

function Recipe_MEG_source_slidingTimeWindow(which_model)

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
toolboxRoot = '/imaging/fj01/latest_toolbox'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code

userOptions = projectOptions();

%%%%%%%%%%%%%%%%%%%%
%% Starting parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%
if userOptions.flush_Queue
    rsa.par.flushQ();
end
if userOptions.run_in_parallel
    p = rsa.par.initialise_CBU_Queue(userOptions);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
userOptions.modelNumber = which_model;
Models = rsa.constructModelRDMs(userOptions);

%%%%%%%%%%%%%%%%%%%
%% Set meta data %%
%%%%%%%%%%%%%%%%%%%
userOptions = rsa.meg.setMetadata_MEG(Models, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sliding time window RoI analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map_type = 'r';

%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Data RDMs %%
%%%%%%%%%%%%%%%%%%%%%%%
tic
rsa.meg.ROI_slidingTimeWindow(userOptions, Models);
toc
%%%%%%%%%%%%%%%%%
%% Permutation %%
%%%%%%%%%%%%%%%%%
tic
if strcmp(userOptions.groupStats,'FFX')
    rsa.meg.FFX_slidingTimeWindow(userOptions,Models);
elseif strcmp(userOptions.groupStats,'RFX')
    rsa.meg.RFX_slidingTimeWindow(userOptions,Models, map_type);
end
toc
%%%%%%%%%%%%%%%%%%%%%
%% Display Results %%
%%%%%%%%%%%%%%%%%%%%%
rsa.meg.showResults_slidingTimeWindow(userOptions, Models, map_type);

%%%%%%%%%%%%%%%%%%%%
%% Sending an email %%
%%%%%%%%%%%%%%%%%%%%
if userOptions.recieveEmail
    rsa.par.setupInternet();
    rsa.par.setupEmail(userOptions.mailto);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stopping parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if userOptions.run_in_parallel;
    % Close the parpool.
    delete(p);
end

%%%%%%%%%%%%%%%%%%%%
%% Delete Selected Directories%%
%%%%%%%%%%%%%%%%%%%%

if (userOptions.deleteTMaps_Dir || userOptions.deleteImageData_Dir || userOptions.deletePerm)
    rsa.util.deleteDir(userOptions, Models);
end


end
