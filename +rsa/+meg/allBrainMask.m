% This function generates dummy index masks fro the case of running whole
% brain analysis (a special mask with all the vertices). Aim of the
% function is to keep homogeneity across all functions and do masking in
% more standard way.
%
% [indexMasks =] allBrainMask(userOptions)
%
% IZ 3/12

function indexMasks = allBrainMask(userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

for chi_i = 1:2
        switch chi_i
            case 1
                chi = 'L';
            case 2
                chi = 'R';
        end%switch:chirality
        
		indexMasks(chi_i).vertices = 1:userOptions.targetResolution; % update IZ 02/12 previously: vertexMask + 1;
		indexMasks(chi_i).timepoints = userOptions.temporalSearchlightLimits;
        indexMasks(chi_i).chirality = chi;
        indexMasls(chi_i).name = ['all_' lower(chi) 'h'];
end

cd(returnHere); % And go back to where you started
