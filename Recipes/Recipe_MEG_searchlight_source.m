% TODO: Documentation
%
% Cai Wingfield 2010-05, 2010-08, 2015-03
% update by Li Su 3-2012, 11-2012
% updated Fawad 12-2013, 02-2014, 10-2014

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
toolboxRoot = 'C:\Users\cai\code\rsagroup-rsatoolbox\'; 
addpath(genpath(toolboxRoot));

userOptions = defineUserOptions();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
models = rsa.constructModelRDMs(userOptions);
% Only using one model at a time for searchlight analysis.
model = models(1);

%%%%%%%%%%%%%%%%%%%%%%
%% Mask preparation %% 
%%%%%%%%%%%%%%%%%%%%%%
usingMasks = ~isempty(userOptions.maskNames);
if usingMasks
    slMasks = rsa.meg.MEGMaskPreparation_source(userOptions);
    % For this searchlight analysis, we combine all masks into one
    slMasks = rsa.meg.combineVertexMasks_source(slMasks, 'combined_mask', userOptions);  
else
    slMasks = rsa.meg.allBrainMask(userOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate adjacency matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSubjects = numel(userOptions.subjectNames);
adjacencyMatrix = rsa.meg.calculateMeshAdjacency(userOptions.targetResolution, userOptions.sourceSearchlightRadius, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Starting parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if userOptions.flush_Queue
    rsa.par.flushQ();
end

if userOptions.run_in_parallel
    p = rsa.par.initialise_CBU_Queue(userOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Searchlight - Brain RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints('Stage 1 - Searchlight Brain RDM Calculation:');
parfor subject_i = 1:nSubjects
    [status, hostname] = system('hostname');
    
    % Work on each hemisphere separately
    for chi = 'LR'
        rsa.util.prints('%s: Working on subject %d, %s side', hostname, subject_i, chi);
    
        % Get subject source data
        [sourceMeshesThisSubjectThisHemi, STCMetadata] = MEGDataPreparation_source( ...
            betaCorrespondence(), ...
            userOptions, ...
            'subject_i', subject_i, ...
            'chi', chi);

        [rMapPaths(subject_i).(chi), RDMPaths(subject_i).(chi)] = rsa.meg.MEGSearchlight_source( ...
            subject_i, ...
            chi, ...
            sourceMeshesThisSubjectThisHemi.(userOptions.subjectNames{subject_i}).(chi), ...
            ...% Use the mask for this hemisphere only
            slMasks([slMasks.chirality] == chi), ...
            model, ...
            adjacencyMatrix, ...
            STCMetadata, ...
            userOptions ...
        );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Random permutation %%
%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(userOptions.groupStats, 'FFX')
    % fixed effect test
    rsa.util.prints('Stage 2 - Fixed Effects Analysis: ');
    tic
    rsa.util.prints('Averaging RDMs across subjects and performing permutation tests to calculate p-values.');
    rsa.meg.FFX_permutation(model, slMasks, userOptions)
    toc
else
    % random effect test
    rsa.util.prints('Stage 2 - Random Effects Analysis: ');
    tic
    rsa.util.prints('Performing permutation tests over all subjects (random effect) to calculate p-values.');
    rsa.meg.RFX_permutation(model,userOptions);
    toc 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sptiotemporal clustering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jobSize = userOptions.jobSize;
number_of_permutations = userOptions.significanceTestPermutations;

tic
rsa.util.prints('Stage 3 - Spatiotemporal Clustering: ' );
parfor j = 1:number_of_permutations/jobSize
    range = (j-1)*jobSize+1:j*jobSize;
    rsa.meg.MEGFindCluster_source(model, range, slMasks, userOptions);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute cluster level p-values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
rsa.util.prints('Stage 4 - Computing cluster level p values: ');
rsa.meg.get_cluster_p_value(model, userOptions);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stopping parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if userOptions.run_in_parallel
    % Close the parpool
    delete(p);
end

%%%%%%%%%%%%%%%%%%%%
%% Delete Selected Directories%%
%%%%%%%%%%%%%%%%%%%%
if (userOptions.deleteTMaps_Dir || userOptions.deleteImageData_Dir || userOptions.deletePerm)
    rsa.util.deleteDir(userOptions, model);
end

%%%%%%%%%%%%%%%%%%%%
%% Sending an email %%
%%%%%%%%%%%%%%%%%%%%
if userOptions.recieveEmail
    rsa.par.setupInternet();
    rsa.par.setupEmail(userOptions.mailto);
end
