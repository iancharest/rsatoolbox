% random_elements = chooseRandom(vector_in[, n_to_select])
%
% Selects random element(s) from the input vector.
%
% EXAMPLE USAGES
%
%     v  = [10 20 30 40 50 60];
%     r1 = chooseRandom(v, 3);  % r1 = [20 10 50]
%     r2 = chooseRandom(v);     % r2 = 40
%
% Cai Wingfield 2015-03
function random_elements = chooseRandom(vector_in, n_to_select)
    %% input checks
    if nargin == 1
        % By default, we select a random element
        n_to_select = 1;
    end

    % a randomised order for the indices of the input vector
    random_order = randperm(numel(vector_in));
    
    % selecting the first indices few gives a random selection
    random_elements = vector_in(random_order(1:n_to_select));
end%function
