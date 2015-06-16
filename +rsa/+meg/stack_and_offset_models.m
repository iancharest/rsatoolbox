% [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data)
%
% Given some dynamic models, and a specified offset in timepoints, produces
% a time-indexed cell array of model stacks suitable for glmfit.
%
% Cai Wingfield 2015-04
function [modelStack, nTimepoints_overlap] = stack_and_offset_models(models, lag_in_timepoints, nTimepoints_data)

    import rsa.*
    import rsa.rdm.*

    [nTimepoints_models, nModels] = size(models);
    
    % We only look at timepoints where the data's timeline overlaps with
    % the models' lag-offset timelines.
    %
    %    (lag)>  |--------------------| lag-offset models
    % |--------------------| data
    %            .         .
    %            |---------|
    %                 ^
    %             only look
    %             in overlap
    nTimepoints_overlap = nTimepoints_data - lag_in_timepoints;
    % There may be a case where we don't have models for the full epoch,
    % even including the overlap, so in this case we ensure that we don't
    % end up looking past the end of the model timeline.
    nTimepoints_overlap = min(nTimepoints_overlap, nTimepoints_models);
    
    model_size = size(models(1,1).RDM);
    
    % Make sure we're using ltv form.
    model_size = size(vectorizeRDM(zeros(model_size)));
    
    % Now at each timepoint we stack the models into a predictor matrix for
    % the GLM.
    % We are only looking in the first bit of the models' timelines, in the
    % places where there is also data.
    for t = 1:nTimepoints_overlap
        for model_i = 1:nModels
            modelStack{t}(model_i, :) = vectorizeRDM(models(t, model_i).RDM);
        end%for:model
    end%for:t
end%function
