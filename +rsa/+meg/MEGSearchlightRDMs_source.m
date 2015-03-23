% MEGSearchlightRDMs_source
%
% Cai Wingfield 2015-03

function [RDMsPath] = MEGSearchlightRDMs_source(subject_i, chi, maskedMeshes, slMask, adjacencyMatrix, STCMetadata, userOptions)

import rsa.*
import rsa.meg.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

subjectName = userOptions.subjectNames{subject_i};

usingMasks = ~isempty(userOptions.maskNames);

%% File paths

RDMsDir = fullfile(userOptions.rootPath, 'RDMs');

if usingMasks
    RDMsFile = ['searchlightRDMs_masked_', subjectName, '-' lower(chi) 'h'];
else
    RDMsFile = ['searchlightRDMs_',        subjectName, '-' lower(chi) 'h'];
end
RDMsPath = fullfile(RDMsDir, RDMsFile);

promptOptions.functionCaller = 'MEGSearchlight_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = RDMsPath;

overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
%% Apply searchlight

if overwriteFlag
    
    prints('Shining RSA searchlight in the source mesh of subject %d of %d...', subject_i, numel(userOptions.subjectNames));
    
    tic;%1
    
    [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions);
    
    searchlightRDMs = rsa.meg.searchlightMappingRDMs_MEG_source(maskedMeshes, slMask, adjacencyMatrix, slSpec, userOptions); %#ok<NASGU>

    %% Saving RDM maps
    
    prints('Saving data RDMs to %s.', RDMsPath);
    gotoDir(RDMsDir);
    save('-v7.3', RDMsPath, 'searchlightRDMs');
    
    %% Done
    t = toc;%1
    prints('That took %s seconds.', t);
else
    prints('Searchlight already applied for subject %s, %s side, skipping it.', subjectName, lower(chi));
end

cd(returnHere); % And go back to where you started

end%function
