% indexMasks_out = combineMasks(indexMasks_in, newMaskName, userOptions)
%
% Combines several masks into a single mask by the union of their vertices.
% Output will have two masks, one called <newMaskName>_l, the other called
% <newMaskName>_r.
%
% Based on code by IZ 2013-03
% Rewritten CW 2015-03
function indexMasks_out = combineVertexMasks_source(indexMasks_in, newMaskName, userOptions)
    % Do left and right hemis separately
    % We'll end up with two masks
    out_mask_i = 1;
    for chi = 'LR'
        
        %% We union the masks together in a vectorised method, for speed
        
        % Masks this hemisphere
        masksThisHemi = indexMasks_in([indexMasks_in.chirality] == chi);
        
        % We want any vertex which is in one of these masks, but without
        % double-counting.
        maskVsThisHemi = unique([masksThisHemi.vertices]);
        
        % Sort the vertices and cap them at the resolution set in
        % userOptions.
        maskVsThisHemi = sort(maskVsThisHemi(maskVsThisHemi <= userOptions.targetResolution));
        
        % Do the same for the timepoints
        maskTsThisHemi = unique([masksThisHemi.timepoints]);
        maskTsThisHemi = sort(maskTsThisHemi);
        
        % Store the unioned bits in a new mask struct
        indexMasks_out(out_mask_i).vertices   = maskVsThisHemi;
        indexMasks_out(out_mask_i).timepoints = maskTsThisHemi;
        indexMasks_out(out_mask_i).chirality  = chi;
        indexMasks_out(out_mask_i).name       = sprintf('%s_%s', newMaskName, lower(chi));
       
        % Make sure we don't overwrite
        out_mask_i = out_mask_i + 1;
    end
end%function
