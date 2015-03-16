%  Recipe_MEG_source
%
% Cai Wingfield 9-2010
% Updated: Isma Zulfiqar 11-2012

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = '/imaging/ls02/toolbox/devel/toolbox'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code
userOptions = defineUserOptions();

models = rsa.constructModelRDMs(userOptions);
model = models(1);

userOptions = rsa.meg.setMetadata_MEG(model, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
nSubjects = numel(userOptions.subjectNames);
for subject = 1:nSubjects
    thisSubject = userOptions.subjectNames{subject};
    fprintf(['Reading MEG source solutions for subject number ' num2str(subject) ' of ' num2str(nSubjects) ': ' thisSubject ':']);    
    sourceMeshes.(thisSubject) = rsa.meg.MEGDataPreparation_source(subject,userOptions.betaCorrespondence, userOptions);
end
indexMasks   = rsa.meg.MEGMaskPreparation_source(userOptions);
maskedMeshes = rsa.meg.MEGDataMasking_source(sourceMeshes, indexMasks, userOptions.betaCorrespondence, userOptions);

%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

RDMs  = rsa.constructRDMs(maskedMeshes, userOptions.betaCorrespondence, userOptions);
RDMs  = rsa.rdm.averageRDMs_subjectSession(RDMs, 'session');
aRDMs = rsa.rdm.averageRDMs_subjectSession(RDMs, 'subject');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rsa.figureRDMs(aRDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1)); % Display the calculated RDMs
rsa.figureRDMs(model, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2)); % Display the models

rsa.MDSConditions(aRDMs, userOptions);
rsa.dendrogramConditions(aRDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

rsa.pairwiseCorrelateRDMs({aRDMs, model}, userOptions);
rsa.MDSRDMs({aRDMs, model}, userOptions);

% Compute distance bar graph comparisons with noise ceiling estimates
userOptions.RDMCorrelationType = 'Kendall_taua';
userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
userOptions.RDMrelatednessThreshold = 0.05;
userOptions.figureIndex = [10 11];
userOptions.RDMrelatednessMultipleTesting = 'FDR';
userOptions.candRDMdifferencesTest = 'subjectRFXsignedRank';
userOptions.candRDMdifferencesThreshold = 0.05;
userOptions.candRDMdifferencesMultipleTesting = 'none';

rsa.compareRefRDM2candRDMs(RDMs(1), model, userOptions);

% fixed effects analysis
rsa.stat.testSignificance({aRDMs}, {model}, userOptions);

% random effects analysis
rsa.stat.testSignificance_RandomEffects(RDMs, model, userOptions)
