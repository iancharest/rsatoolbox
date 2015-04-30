% TODO: Documentation
%
% Cai Wingfield 2010-05, 2010-08, 2015-03, 2015-04
% update by Li Su 3-2012, 11-2012
% updated Fawad 12-2013, 02-2014, 10-2014

toolboxRoot = '/imaging/cw04/code/rsagroup-rsatoolbox/';
addpath(genpath(toolboxRoot));

userOptions = defineUserOptions();

rsa.util.prints('Starting RSA analysis "%s".', userOptions.analysisName);


%% %%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Preparing model RDMs...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%

models = rsa.constructModelRDMs(userOptions);


%% %%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Preparing masks...');
%%%%%%%%%%%%%%%%%%%%%%

usingMasks = ~isempty(userOptions.maskNames);
if usingMasks
    slMasks = rsa.meg.MEGMaskPreparation_source(userOptions);
    % For this searchlight analysis, we combine all masks into one
    slMasks = rsa.meg.combineVertexMasks_source(slMasks, 'combined_mask', userOptions);  
else
    slMasks = rsa.meg.allBrainMask(userOptions);
end


%% Compute some constats
nSubjects = numel(userOptions.subjectNames);
adjacencyMatrices = rsa.meg.calculateMeshAdjacency(userOptions.targetResolution, userOptions.sourceSearchlightRadius, userOptions, 'hemis', 'LR');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Starting parallel toolbox...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if userOptions.flush_Queue
    rsa.par.flushQ();
end

if userOptions.run_in_parallel
    p = rsa.par.initialise_CBU_Queue(userOptions);
end


%% %%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Loading brain data...');
%%%%%%%%%%%%%%%%%%%%%

[meshPaths, STCMetadatas] = rsa.meg.MEGDataPreparation_source( ...
    betaCorrespondence(), ...
    userOptions, ...
    'mask', slMasks);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Searchlight Brain RDM Calculation...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RDMsPaths, slSTCMetadatas] = rsa.meg.MEGSearchlightRDMs_source( ...
    meshPaths, ...
    slMasks, ...
    ...% Assume that both hemis' adjacency matrices are the same so only use one.
    adjacencyMatrices.L, ...
    STCMetadatas, ...
    userOptions);


%% %%%%%
rsa.util.prints( ...
    'Averaging searchlight RDMs...');
%%%%%%%%

averageRDMPaths = rsa.meg.averageSearchlightRDMs(RDMsPaths, userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'GLM-fitting models to searchlight RDMs...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for subject_i = 1:numel(userOptions.subjectNames)
    subjectName = userOptions.subjectNames{subject_i};
    [glm_paths.(subjectName), lagSTCMetadatas.(subjectName)] = rsa.meg.searchlight_dynamicGLM_source( ...
        RDMPaths(subject_i), ...
        models, ...
        slSTCMetadatas, ...
        userOptions, ...
        'lag', 100, ...
        'file-prefix', [subjectName '-']);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Thresholding GLM values...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p_paths, p_median_paths] = rsa.meg.searchlight_GLM_permutation_source( ...
    averageRDMPaths, ...
    glm_paths, ...
    models, ...
    slSTCMetadatas, ...
    lagSTCMetadatas, ...
    30, ...
    userOptions);

[thresholded_glm_paths] = rsa.meg.searchlight_GLM_threshold_source( ...
    glm_paths, ...
    p_paths, ...
    p_median_paths, ...
    lagSTCMetadatas, ...
    0.05);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Cleaning up...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close the parpool
if userOptions.run_in_parallel
    delete(p);
end

% Sending an email
if userOptions.recieveEmail
    rsa.par.setupInternet();
    rsa.par.setupEmail(userOptions.mailto);
end

rsa.util.prints( ...
    'RSA COMPLETE!');
