% MEGSearchlight_source
%
% MEGSearchlight_source(subjectNumber, sourceMeshesThisSubject, indexMasks, Models, userOptions)
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

% TODO: documentation

function MEGSearchlight_source(subjectNumber, chi, sourceMeshesThisSubjectThisHemi, indexMask, Models, adjacencyMatrix, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

subject = userOptions.subjectNames{subjectNumber};
nSubjects = numel(userOptions.subjectNames);

% TODO: this is a weird thing to have in userOptions.
modelNumber = userOptions.modelNumber;
modelName = spacesToUnderscores(Models(modelNumber).name);

if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

if userOptions.maskingFlag
    MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, '_', subject '_masked'];
else
    MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, '_', subject];
end

promptOptions.functionCaller = 'MEGSearchlight_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-', lower(chi), 'h.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    
    prints('Shining RSA searchlights...');
    gotoDir(fullfile(userOptions.rootPath, 'Maps'), modelName);
    
    tic;%1
    
    prints(['\tSearching in the source meshes of subject ' num2str(subjectNumber) ' of ' num2str(nSubjects) ':']);
    
    %% Mask timepoints
    if userOptions.nSessions == 1
        maskedMesh = sourceMeshesThisSubjectThisHemi(:, timeIndices(1):timeIndices(2), :); % (vertices, timePointes, conditions)
    else
        maskedMesh = sourceMeshesThisSubjectThisHemi(:, timeIndices(1):timeIndices(2), :, :); % (vertices, timePointes, conditions, sessions)
    end

    %% Apply searchlight
    if userOptions.partial_correlation
        % It says searchlightRDMs are unused, but actually they're saved.
        [thisSubjectRs, thisSubjectPs, searchlightRDMs] = searchlightMapping_MEG_source(maskedMesh, indexMask, Models(modelNumber), Models([userOptions.partial_modelNumber{:}]), adjacencyMatrix, userOptions); %#ok<ASGLU>
    else
        [thisSubjectRs, thisSubjectPs, searchlightRDMs] = searchlightMapping_MEG_source(maskedMesh, indexMask, Models(modelNumber), [], adjacencyMatrix, userOptions); %#ok<ASGLU>
    end

    rMetadataStruct = userOptions.STCmetaData;

    rMetadataStruct.data = thisSubjectRs(:,:,modelNumber);

    %% Saving r-maps and p-maps
    if userOptions.maskingFlag
        outputRFilename = [fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_rMesh_' modelName '_' subject ]) '_masked' '-' lower(chi) 'h.stc'];
    else
        outputRFilename = [fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_rMesh_' modelName '_' subject]) '-' lower(chi) 'h.stc'];
    end
    
    % Write the r and p maps
    mne_write_stc_file1(outputRFilename, rMetadataStruct);

    %% Saving the searchlight RDMs
    prints('Saving data RDMs for combined mask: ');
    filepath = 'searchlightRDMs_';
    if userOptions.maskingFlag
        filepath = [filepath 'masked_'];
    end
    gotoDir(userOptions.rootPath, 'RDMs');
    save('-v7.3', [filepath, userOptions.subjectNames{subjectNumber}, '-', lower(chi), 'h'], 'searchlightRDMs');

    % Free up some memory
    clear thisSubjectRs thisSubjectPs pMetadataStruct searchlightRDMs maskedMesh;

    t = toc;%1
    prints('That took %s seconds.', t);
else
    prints('Searchlight already applied, skipping it.');
end

cd(returnHere); % And go back to where you started

end%function
