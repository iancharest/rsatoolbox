% [modelStack] = stack_and_trim_models(models, nTimepoints_overlap)
%
% Given some dynamic models, and a specified offset in timepoints, produces
% a time-indexed cell array of model stacks suitable for glmfit, trims
% everything which won't overlap with the data timeline.
%
% Cai Wingfield 2015-04, 2015-08
function [modelStack] = stack_and_trim_models(models, nTimepoints_overlap)

    import rsa.*
    import rsa.rdm.*

    [nTimepoints_models, nModels] = size(models);
    
    % We only look at timepoints where the data's timeline overlaps with
    % the models' lag-offset timelines.
    %
    %    (lag)>  |----------(---trim---)| lag-offset models
    % |--------------------| data
    %            .         .
    %            |---------|
    %                 ^
    %             only look
    %             in overlap
    
    model_size = size(models(1,1).RDM);
    
    % Make sure we're using ltv form.
    % (value used for debugging only rn)
    model_size = size(vectorizeRDM(zeros(model_size)));
    
    % Now at each timepoint we stack the models into a predictor matrix for
    % the GLM.
    % We are only looking in the first bit of the models' timelines, in the
    % places where there is also data.
    for t = 1:nTimepoints_overlap
        for model_i = 1:nModels
            modelStack{t}(model_i, :) = vectorizeRDM(models( ...
                t, ...
                model_i).RDM);
        end%for:model
    end%for:t
end%function
