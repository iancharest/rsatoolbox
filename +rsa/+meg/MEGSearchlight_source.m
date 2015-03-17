% MEGSearchlight_source
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

% TODO: Documentation
% TODO: Partial model numbers vs model number - make this coherent

function MEGSearchlight_source(subject_i, chi, sourceMeshes, indexMask, model, partialModels, adjacencyMatrix, STCMetadata, userOptions)

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

modelName = spacesToUnderscores(model.name);

if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

if usingMasks
    MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, '_', subjectName '_masked'];
else
    MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, '_', subjectName];
end

promptOptions.functionCaller = 'MEGSearchlight_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-', lower(chi), 'h.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    
    prints('Shining RSA searchlight in the source mesh of subject %d of %d...', subject_i, nSubjects);
    
    prints('Shining RSA searchlights...');
    gotoDir(fullfile(userOptions.rootPath, 'Maps', modelName));
    
    tic;%1
    
    [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions);
    
    %% Apply searchlight
    if userOptions.partial_correlation
        [thisSubjectRs, searchlightRDMs] = searchlightMapping_MEG_source(sourceMeshes, indexMask, model, partialModels, adjacencyMatrix, slSpec, userOptions); %#ok<ASGLU>
    else
        [thisSubjectRs, searchlightRDMs] = searchlightMapping_MEG_source(sourceMeshes, indexMask, model, [], adjacencyMatrix, slSpec, userOptions); %#ok<ASGLU>
    end
    clear sourceMeshesThisSubjectThisHemi;

    rSTCStruct = slSTCMetadata;
    rSTCStruct.data = thisSubjectRs(:,:);

    %% Saving r-maps
    if usingMasks
        outputRFilename = [fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_rMesh_' modelName '_' subjectName ]) '_masked' '-' lower(chi) 'h.stc'];
    else
        outputRFilename = [fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_rMesh_' modelName '_' subjectName]) '-' lower(chi) 'h.stc'];
    end
    
    % Write the r-maps
    prints('Writing r-map %s.', outputRFilename);
    mne_write_stc_file1(outputRFilename, rSTCStruct);

    %% Saving the searchlight RDMs
    
    if usingMasks
        filepath = 'searchlightRDMs_masked_';
    else
        filepath = 'searchlightRDMs_';
    end
    gotoDir(userOptions.rootPath, 'RDMs');
    saveLocation = [filepath, userOptions.subjectNames{subject_i}, '-', lower(chi), 'h'];
    prints('Saving data RDMs for combined mask to %s.', saveLocation);
    save('-v7.3', saveLocation, 'searchlightRDMs');

    t = toc;%1
    prints('That took %s seconds.', t);
else
    prints('Searchlight already applied, skipping it.');
end

cd(returnHere); % And go back to where you started

end%function
