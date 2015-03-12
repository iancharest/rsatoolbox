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
        maskVertices = [];
        for mask_i = 1:numel(indexMasks_in)
            if strcmpi(indexMasks_in(mask_i).chirality, chi)
                maskVertices = union(maskVertices, indexMasks(mask_i).vertices);
            end
        end
        maskVertices = sort(maskVertices(maskVertices <= userOptions.nVertices));
        % TODO: this isn't coming from the individual masks... is that ok?
        timeIndices = userOptions.dataPointsSearchlightLimits;
        
        % Store the unioned bits in a new mask struct
        indexMasks_out(out_mask_i).vertices   = maskVertices;
        indexMasks_out(out_mask_i).timepoints = timeIndices;
        indexMasks_out(out_mask_i).chirality  = chi;
        indexMasks_out(out_mask_i).name       = sprintf('%s_%s', newMaskName, lower(chi));
       
        % Make sure we don't overwrite
        out_mask_i = out_mask_i + 1;
    end
end%function
