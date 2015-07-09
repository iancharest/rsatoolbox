% Given a binary vector representing the presence of absence of a
% particular feature for each condition, this function returns a model RDM
% based on that condition.
%
% Returns the RDM in ltv-form.
%
% CW 2015-05
function rdm = binary_categorical_rdm(v)

    n_conditions = numel(v);
    
    rdm = zeros(n_conditions, n_conditions);
    for condition_1 = 1:n_conditions-1
        for condition_2 = condition_1+1:n_conditions
            rdm(condition_2, condition_1) = (v(condition_1) ~= v(condition_2));
        end
    end
    
    rdm = squareform(rdm);

end%function
