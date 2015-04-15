% searchlightMapping_source
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

% TODO: Documentation
% TODO: Partial model numbers vs model number - make this coherent

function [mapsPath] = searchlightMapping_source(subject_i, chi, RDMPath, slMask, model, partialModels, adjacencyMatrix, STCMetadatas, userOptions)

import rsa.*
import rsa.meg.*
import rsa.rdm.*
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

mapsDir = fullfile(userOptions.rootPath, 'Maps', modelName);
if usingMasks
    mapsFile = [userOptions.analysisName '_rMesh_' modelName '_' subjectName '_masked-' lower(chi) 'h.stc'];
else
    mapsFile = [userOptions.analysisName '_rMesh_' modelName '_' subjectName '-' lower(chi) 'h.stc'];
end
mapsPath = fullfile(mapsDir, mapsFile);

promptOptions.functionCaller = 'searchlightMapping_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = mapsPath;

overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
%% Apply searchlight

if overwriteFlag
    
    prints('Shining RSA searchlight in the %s source mesh of subject %d of %d...', lower(chi),  subject_i, nSubjects);
    
    gotoDir(fullfile(userOptions.rootPath, 'Maps', modelName));
    
    tic;%1
    
    [slSpecs, slSTCMetadatas] = getSearchlightSpec(STCMetadatas, userOptions);
    
    searchlightRDMs = directLoad(RDMPath, 'searchlightRDMs');
    
    if userOptions.partial_correlation
        thisSubjectRs = searchlightMapping_MEG_source(searchlightRDMs, slMask, model, partialModels, adjacencyMatrix, slSpecs.(chi), userOptions); %#ok<ASGLU>
    else
        thisSubjectRs = searchlightMapping_MEG_source(searchlightRDMs, slMask, model, [], adjacencyMatrix, slSpecs.(chi), userOptions); %#ok<ASGLU>
    end
    clear sourceMeshesThisSubjectThisHemi;

    rSTCStruct          = slSTCMetadatas.(chi);
    rSTCStruct.vertices = slMask.vertices;
    % thisSubjectRs contains only the data inside the mask, but since the
    % vertices are stored in this struct, that should be ok.
    rSTCStruct.data     = thisSubjectRs(:,:);

    %% Saving r-maps and RDM maps
    
    prints('Writing r-map %s.', mapsPath);
    gotoDir(mapsDir);
    mne_write_stc_file1(mapsPath, rSTCStruct);
    
    %% Done
    t = toc;%1
    prints('That took %s seconds.', t);
else
    prints('Searchlight already applied, skipping it.');
end

cd(returnHere); % And go back to where you started

end%function

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

% [smm_rs, searchlightRDMs] = searchlightMapping_MEG_source(singleSubjectMesh, indexMask, modelRDM, partialModelRDMs, adjacencyMatrix, slSpec, userOptions)
%
% Based on Li Su's script
% CW 2010-05, 2015-03
% updated by Li Su 3-2012

function [smm_rs, searchlightRDMs] = searchlightMapping_MEG_source(searchlightRDMs, indexMask, modelRDM, partialModelRDMs, adjacencyMatrix, slSpecs, userOptions)

import rsa.*
import rsa.meg.*
import rsa.rdm.*
import rsa.stat.*
import rsa.util.*

modelRDM_utv = squeeze(unwrapRDMs(vectorizeRDMs(modelRDM)));

if userOptions.partial_correlation
    % TODO: should be transposed?
    control_for_modelRDMs = unwrapRDMs(vectoriseRDMs(partialModelRDMs));
end

[nVertices_masked, nTimePoints_masked] = size(searchlightRDMs);

% The number of positions the sliding window will take.
nWindowPositions = size(slSpecs.windowPositions, 1);

%% map the volume with the searchlight

% Preallocate looped matrices for speed
smm_rs = zeros(nVertices_masked, nWindowPositions);

% For display purposes
nVertsSearched = 0;

% Search the vertices
for v_i = 1:numel(indexMask.vertices)
    % v_i loops through the *indices of* vertices in the mask
    % v is the vertex itself
    
    v = indexMask.vertices(v_i);
    
    % Search through time
    window_i = 0;
    for window = slSpecs.windowPositions'
        % thisWindow is the indices of timepoints in each window
        thisWindow = window(1):window(2);
        window_i = window_i + 1;
        
        searchlightRDM = searchlightRDMs(v, window).RDM;
        
        % TODO: Refactor this into general method so it can be used
        % TODO: anywhere (this is being done on another branch)
        if strcmpi(userOptions.RDMCorrelationType, 'Kendall_taua')
            rs = rankCorr_Kendall_taua(searchlightRDM', modelRDM_utv');
        elseif userOptions.partial_correlation
            % TODO: Consider partialcorr with Kendall's tau
            rs = partialcorr(searchlightRDM', modelRDM_utv', control_for_modelRDMs', 'type', userOptions.RDMCorrelationType, 'rows','pairwise');
        else
            rs = corr(searchlightRDM', modelRDM_utv', 'type', userOptions.RDMCorrelationType, 'rows', 'pairwise');
        end
        
        % Store results to be retured.
        smm_rs(v_i, window_i) = rs;
        
    end%for:window
    
    % Indicate progress every once in a while...
    nVertsSearched = nVertsSearched + 1;
    if mod(nVertsSearched, 200) == 0
        prints('%d vertices searched: %d%% complete', nVertsSearched, floor((nVertsSearched / numel(indexMask.vertices)) * 100));
    end%if
    
end%for:v

if userOptions.fisher
    smm_rs = fisherTransform(smm_rs);
end%if

end%function

