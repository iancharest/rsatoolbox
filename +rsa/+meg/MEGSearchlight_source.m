% MEGSearchlight_source
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

% TODO: Documentation
% TODO: Partial model numbers vs model number - make this coherent

function [mapsPath, RDMsPath] = MEGSearchlight_source(subject_i, chi, maskedMeshes, slMask, model, partialModels, adjacencyMatrix, STCMetadata, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

subjectName = userOptions.subjectNames{subject_i};
nSubjects = numel(userOptions.subjectNames);

usingMasks = ~isempty(userOptions.maskNames);

if userOptions.partial_correlation
    modelName = [spacesToUnderscores(model.name), '_partialCorr'];
else
    modelName = spacesToUnderscores(model.name);
end

%% File paths

RDMsDir = fullfile(userOptions.rootPath, 'RDMs');
mapsDir = fullfile(userOptions.rootPath, 'Maps', modelName);
if usingMasks
    RDMsFile = ['searchlightRDMs_masked_', userOptions.subjectNames{subject_i}, '-' lower(chi) 'h'];
    mapsFile = [userOptions.analysisName '_rMesh_' modelName '_' subjectName '_masked-' lower(chi) 'h.stc'];
else
    RDMsFile = ['searchlightRDMs_', userOptions.subjectNames{subject_i}, '-' lower(chi) 'h'];
    mapsFile = [userOptions.analysisName '_rMesh_' modelName '_' subjectName '-' lower(chi) 'h.stc'];
end
RDMsPath = fullfile(RDMsDir, RDMsFile);
mapsPath = fullfile(mapsDir, mapsFile);

promptOptions.functionCaller = 'MEGSearchlight_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = RDMsPath;
promptOptions.checkFiles(2).address = mapsPath;

overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
%% Apply searchlight

if overwriteFlag
    
    prints('Shining RSA searchlight in the source mesh of subject %d of %d...', subject_i, nSubjects);
    
    prints('Shining RSA searchlights...');
    gotoDir(fullfile(userOptions.rootPath, 'Maps', modelName));
    
    tic;%1
    
    [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions);
    
    if userOptions.partial_correlation
        [thisSubjectRs, searchlightRDMs] = searchlightMapping_MEG_source(maskedMeshes, slMask, model, partialModels, adjacencyMatrix, slSpec, userOptions); %#ok<ASGLU>
    else
        [thisSubjectRs, searchlightRDMs] = searchlightMapping_MEG_source(maskedMeshes, slMask, model, [], adjacencyMatrix, slSpec, userOptions); %#ok<ASGLU>
    end
    clear sourceMeshesThisSubjectThisHemi;

    rSTCStruct          = slSTCMetadata;
    rSTCStruct.vertices = slMask.vertices;
    % thisSubjectRs contains only the data inside the mask, but since the
    % vertices are stored in this struct, that should be ok.
    rSTCStruct.data     = thisSubjectRs(:,:);

    %% Saving r-maps and RDM maps
    
    prints('Writing r-map %s.', mapsPath);
    gotoDir(mapsDir);
    mne_write_stc_file1(mapsPath, rSTCStruct);

    prints('Saving data RDMs for combined mask to %s.', RDMsPath);
    gotoDir(RDMsDir);
    save('-v7.3', RDMsPath, 'searchlightRDMs');
    
    %% Done
    t = toc;%1
    prints('That took %s seconds.', t);
else
    prints('Searchlight already applied, skipping it.');
end

cd(returnHere); % And go back to where you started

end%function
