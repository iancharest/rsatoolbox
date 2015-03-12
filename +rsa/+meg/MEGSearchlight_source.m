% MEGSearchlight_source
%
% MEGSearchlight_source(subjectNumber, sourceMeshesThisSubject, indexMasks, Models, userOptions)
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

function MEGSearchlight_source(subjectNumber, sourceMeshesThisSubject, indexMasks, Models, adjacencyMatrix, userOptions)

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
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-lh.stc']);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-rh.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    
    prints('Shining RSA searchlights...');
    gotoDir(userOptions.rootPath, 'Maps');
    gotoDir(fullfile(userOptions.rootPath, 'Maps'), modelName);
    
    tic;%1
    
    prints(['\tSearching in the source meshes of subject ' num2str(subjectNumber) ' of ' num2str(nSubjects) ':']);
    
    for chirality = 1:2
        switch chirality
            case 1
                chi = 'L';
            case 2
                chi = 'R';
        end%switch:chirality
        
        %% Mask timepoints
        if userOptions.nSessions == 1
            maskedMesh = sourceMeshesThisSubject.(chi)(:, timeIndices(1):timeIndices(2), :); % (vertices, timePointes, conditions)
        else
            maskedMesh = sourceMeshesThisSubject.(chi)(:, timeIndices(1):timeIndices(2), :, :); % (vertices, timePointes, conditions, sessions)
        end
        
        %% Apply searchlight
        if userOptions.partial_correlation
            [thisSubjectRs.(chi), thisSubjectPs.(chi), searchlightRDMs] = searchlightMapping_MEG_source(maskedMesh, Models(modelNumber), Models([userOptions.partial_modelNumber{:}]), adjacencyMatrix, userOptions);
        else
            [thisSubjectRs.(chi), thisSubjectPs.(chi), searchlightRDMs] = searchlightMapping_MEG_source(maskedMesh, Models(modelNumber), [], adjacencyMatrix, userOptions);
        end
        
        rMetadataStruct = userOptions.STCmetaData;
        pMetadataStruct = userOptions.STCmetaData;
        
        rMetadataStruct.data = thisSubjectRs.(chi)(:,:,modelNumber);
        pMetadataStruct.data = thisSubjectPs.(chi)(:,:,modelNumber);
        
        %% Saving r-maps and p-maps
        outputRFilename = fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_rMesh_' modelName '_' subject ]);
        outputPFilename = fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_pMesh_' modelName '_' subject ]);
        if userOptions.maskingFlag
            outputRFilename = [outputRFilename '_masked'];
            outputPFilename = [outputPFilename '_masked'];
        end
        mne_write_stc_file1([outputRFilename '-' lower(chi) 'h.stc'], rMetadataStruct);
        mne_write_stc_file1([outputPFilename '-' lower(chi) 'h.stc'], pMetadataStruct);
        
        %% Saving the searchlight RDMs
        prints('Saving data RDMs for combined mask: ');
        filepath = 'searchlightRDMs_';
        if userOptions.maskingFlag
            filepath = [filepath 'masked_']; %#ok<AGROW>
        end
        gotoDir(userOptions.rootPath, 'RDMs');
        save('-v7.3', [filepath, userOptions.subjectNames{subjectNumber},'-',lower(chi),'h'], 'searchlightRDMs');
        
        userOptions = rmfield(userOptions, 'maskIndices');
        userOptions = rmfield(userOptions, 'chi');
        clear thisSubjectRs thisSubjectPs pMetadataStruct searchlightRDMs sourceMeshes.(chi) maskedMesh;
        
    end%for:chirality
    
    % Print the elapsed time for this subject
    if chirality == 1
        fprintf('\n\t\t\t\t\t\t\t\t');
    else
        t = toc;%1
        prints([': [' num2str(ceil(t)) 's]']);
    end%if
    
else
    prints('Searchlight already applied, skip....\n');
end

cd(returnHere); % And go back to where you started

end%function
