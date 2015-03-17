% Recipe_MEG_searchlight_sensor
%
% Cai Wingfield 3-2010, 8-2010

function Recipe_MEG_searchlight_sensor(which_model)

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
toolboxRoot = '/imaging/iz01/test/toolbox/devel/toolbox/'; addpath(genpath(toolboxRoot));
userOptions = projectOptions();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = rsa.constructModelRDMs(userOptions);
model = model(which_model);

%%%%%%%%%%%%%%%%%%%%
%% Set Meta-data %%
%%%%%%%%%%%%%%%%%%%%
% TODO: This isn't done anymore. Make sure that it's excised.
userOptions = rsa.meg.setMetadata_MEG(model, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
sensorImages = rsa.meg.MEGDataPreparation_sensor(userOptions);

%%%%%%%%%%%%%%%%%
%% Searchlight %%
%%%%%%%%%%%%%%%%%
rsa.meg.MEGSearchlight_sensor(sensorImages, model, userOptions);

end
