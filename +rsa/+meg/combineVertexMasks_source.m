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
        % We want the unique vertices
        maskVertices = unique( ...
            ...% which are present in those input masks...
            [indexMasks_in( ...
                ...% where...
                find( ...
                    ...% the chirality is the same as the hemisphere we're looking at
                    [indexMasks_in.chirality] == chi)).vertices] ...
        );
        
        % Sort the vertices and cap them at the resolution set in
        % userOptions.
        maskVertices = sort(maskVertices(maskVertices <= userOptions.targetResolution));
        
        % Set the timeIndices also
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
