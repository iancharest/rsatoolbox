% searchlightAdjacencies = calculateMeshAdjacency(nVertices, searchlightRadius_mm, userOptions, ['hemi', 'L'|'R'|'LR'])
%
% All credit to Su Li and Andy Thwaites for working out how to do this and writing the original implementation
% CW 5-2010, last updated by Li Su - 1 Feb 2012
function searchlightAdjacencies = calculateMeshAdjacency(nVertices, searchlightRadius_mm, userOptions, varargin)

    import rsa.*
    import rsa.fig.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.sim.*
    import rsa.spm.*
    import rsa.stat.*
    import rsa.util.*

    %% Parse inputs

    % 'hemis'
    nameHemis = 'hemis';
    validHemis = {'L', 'l', 'R', 'r', 'LR', 'lr', 'RL', 'rl'};
    checkHemis = @(x) (any(validatestring(x, validHemis)));
    defaultHemis = 'L';

    % Set up parser
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.StructExpand = false;

    % Parameters
    addParameter(ip, nameHemis, defaultHemis, checkHemis);

    % Parse the inputs
    parse(ip, {varargin{:});

    % Which hemisphere are we using?
    % Use upper to keep the convention consistent and easy to read.
    chis = upper(ip.Results.hemis);

    returnHere = pwd;


    %% Set up some constants

    % The maximum possible per-hemisphere resolution.
    MAX_VERTICES = 40968;

    % Freesurfer resolution.
    FS_RESOLUTION = 1.25; % mm between freesurfer vertices

    % Work out some derived values which don't change per hemisphere

    downsampleRate = log(ceil(MAX_VERTICES / nVertices)) / log(4) * 4; 
    searchlightRadius_freesurfer = searchlightRadius_mm / FS_RESOLUTION;
    searchlightCircleRadii_MNE = ceil(downsampleRate * (1:(searchlightRadius_freesurfer / downsampleRate)));
    
    gotoDir(userOptions.rootPath, 'ImageData');

    for chi = chis

        matrixFilename = sprintf('%s_vertexAdjacencyTable_radius-%dmm_%d-verts-%sh.mat', userOptions.analysisname, searchlightRadius_mm, nVertices, lower(chi));

        if ~exist(matrixFilename, 'file')

            prints('The file "%s" doesn''t exist yet, creating it...', matrixFilename);

            % Building a hash table to store adjacency information of all vertexs. 
            % The resulting hash table is ht_*
            hashTable = findAdjacentVerts(userOptions.averageSurfaceFiles.(chi));

            prints('Building vertex adjacency matrix...');

            parfor currentSearchlightCentre = 1:nVertices

                % Print feedback every once in a while.
                if mod(currentSearchlightCentre, floor(nVertices/11)) == 0
                    prints('Working on vertex %d of %d', currentSearchlightCentre, nVertices, floor(100*(currentSearchlightCentre/nVertices)));
                end

                verticesWithinSearchlight = [];

                for rMNE = 1:numel(searchlightCircleRadii_MNE)
                    freesurferVerticesWithinThisMNERadius = getadjacent(num2str(currentSearchlightCentre),searchlightCircleRadii_MNE(rMNE),hashTable);
                    verticesWithinSearchlight = [verticesWithinSearchlight; freesurferVerticesWithinThisMNERadius(freesurferVerticesWithinThisMNERadius <= nVertices)]; % By removing any which are greater than nVertices, we effectively downsample by the necessary ammount.  This seems a little too clever to work? < or <=?
                end

                searchlightAdjacency(currentSearchlightCentre,1:numel(verticesWithinSearchlight)) = verticesWithinSearchlight';
            end

            prints('Done!');

            % Force nans where zeros are.
            searchlightAdjacency(searchlightAdjacency == 0) = NaN;

            % Save this matrix
            save(matrixFilename, 'searchlightAdjacency');

        else

            prints('The file "%s" has already been created, loading it...', matrixFilename);

            searchlightAdjacency = directLoad(matrixFilename, 'searchlightAdjacency');

        end
        
        searchlightAdjacencies.(chi) = searchlightAdjacency;
        clear searchlightAdjacency;
        
    end%for:chis

    cd(returnHere);

end%function

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

% this function returns a hash table containing the adjacent vertexes for
% each vertex in the brain based on freesurfer cortical mash.
% 
% created by Li Su and Andy Thwaites, last updated by Li Su 01 Feb 2012
function ht = findAdjacentVerts(path)

    import rsa.*
    import rsa.fig.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.sim.*
    import rsa.spm.*
    import rsa.stat.*
    import rsa.util.*

    % addpath /opt/mne/matlab/toolbox/ % CW: path doesn't exist.

    [verts,faces] = mne_read_surface(path);
    numberOfVerts = max(faces);
    numberOfFaces = size(faces);

    ht = java.util.Hashtable;

    facesduplicate = zeros(length(faces)*3, 3);

    for i = 1:length(faces)
        q = length(faces);
        % disp(num2str(i));
        facesduplicate(i,1:3) = [faces(i,1) faces(i,2) faces(i,3)];
        facesduplicate(i+q,1:3) = [faces(i,2) faces(i,1) faces(i,3)];
        facesduplicate(i+(q*2),1:3) = [faces(i,3) faces(i,2) faces(i,1)];
    end

    sortedfaces = sortrows(facesduplicate,1);

    thisface = 1;
    adjacent = [];
    for i = 1:length(sortedfaces)
        % disp(num2str(i));
        face = sortedfaces(i,1);
        if  (face == thisface)
            key = num2str(face);
            adjacent = [adjacent sortedfaces(i,2)];
            adjacent = [adjacent sortedfaces(i,3)];
        else
            unad = unique(adjacent);
            ht.put(key,unad);
            adjacent = [];
            thisface = face;

            % now continue as normall
            key = num2str(face);
            adjacent = [adjacent sortedfaces(i,2)];
            adjacent = [adjacent sortedfaces(i,3)];
        end
    end

end%function

% by Li Su and Andy Thwaites
function [adjacents, passed] = getadjacent(str1, int, hashtab)

    import rsa.*
    import rsa.fig.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.sim.*
    import rsa.spm.*
    import rsa.stat.*
    import rsa.util.*

    adjacentsbelow = [];
    adjacents = [];
    passed = [];
    
    if int==1
       adjacents = hashtab.get(str1);
       passed = [1];
    else
       [adjacentsbelow, passed] = getadjacent(str1, int-1, hashtab);
       for j = 1:length(adjacentsbelow)
          adjacents = [adjacents; hashtab.get(num2str(adjacentsbelow(j)))];
       end
       adjacents = unique(adjacents);
       passed = [passed; adjacentsbelow];
       for j = length(adjacents):-1:1
          if(any(find(passed == adjacents(j))))
              adjacents(j)=[];
          end
       end
    end

    adjacents = unique(adjacents);
    passed = unique(passed);

end%function

