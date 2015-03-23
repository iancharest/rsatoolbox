function searchlightAdjacency = calculateMeshAdjacency(nVertices, searchlightRadius_mm, userOptions)

% searchlightAdjacency = calculateMeshAdjacency(nVertices, searchlightRadius_mm, userOptions)
%
% All credit to Su Li and Andy Thwaites for working out how to do this and writing the original implementation
% CW 5-2010, last updated by Li Su - 1 Feb 2012

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd;

MAX_VERTICES = 40968;

downsampleRate = log(ceil(MAX_VERTICES / nVertices)) / log(4) * 4; 
freesurferResolution = 1.25; % mm between freesurfer vertices

searchlightRadius_freesurfer = searchlightRadius_mm / freesurferResolution;

searchlightCircleRadii_MNE = ceil(downsampleRate*(1:(searchlightRadius_freesurfer/downsampleRate)));

matrixFilename = [userOptions.analysisName '_vertexAdjacencyTable_radius-' num2str(searchlightRadius_mm) 'mm_' num2str(nVertices) '-vertices.mat'];

gotoDir(userOptions.rootPath, 'ImageData');

if ~exist(matrixFilename, 'file')

	fprintf('The file "%s" doesn''t exist yet, creating it...', matrixFilename);
	
	% building a hash table to store adjacency information of all vertexs. 
	hashTableL = findAdjacentVerts(userOptions.averageSurfaceFile); % the resulting hash table is ht_*
	% I might be wrong, but I have discovered that the adjacent
	% vertex indexing is the same across hemispheres, so I only do
	% it for the left. Li Su
	
	prints('Building vertex adjacency matrix...');
		
	parfor currentSearchlightCentre = 1:nVertices
	
		if mod(currentSearchlightCentre,floor(nVertices/11)) == 0
			prints(['   Working on the vertex ' num2str(currentSearchlightCentre) ' of ' num2str(nVertices) ': ' num2str(floor(100*(currentSearchlightCentre/nVertices))) '%%']);
		end
		
		verticesWithinSearchlight = [];
			
		for rMNE = 1:numel(searchlightCircleRadii_MNE)
			freesurferVerticesWithinThisMNERadius = getadjacent(num2str(currentSearchlightCentre),searchlightCircleRadii_MNE(rMNE),hashTableL);
			verticesWithinSearchlight = [verticesWithinSearchlight; freesurferVerticesWithinThisMNERadius(freesurferVerticesWithinThisMNERadius <= nVertices)]; % By removing any which are greater than nVertices, we effectively downsample by the necessary ammount.  This seems a little too clever to work? < or <=?
		end
		
		searchlightAdjacency(currentSearchlightCentre,1:numel(verticesWithinSearchlight)) = verticesWithinSearchlight';
	end
	
	prints('      Done!');
	
	searchlightAdjacency(searchlightAdjacency == 0) = NaN;
	
	% Save this matrix
	
	cd(fullfile(userOptions.rootPath, 'ImageData'));
	
	save(matrixFilename, 'searchlightAdjacency');

else

	prints('The file "%s" has already been created, loading it...', matrixFilename);

	load(matrixFilename);

end

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

