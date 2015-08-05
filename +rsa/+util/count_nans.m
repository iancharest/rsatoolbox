% n_nans = count_nans(X)
%
% Counts the number of nans in a multidimensional array.
%
% CW 2015-08
function n_nans = count_nans(X)
    n_nans = sum(isnan(X(:)));
end%function
