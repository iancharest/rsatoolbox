% element_list = elements(matrix_in)
%
% Lists the elements in a n-dimensional array.
%
% EG:
% >> x = randn(2,2,2,2);
% >> elements(a)
%  ans =
%   -0.2099
%   -1.6989
%    0.6076
%   -0.1178
%    0.6992
%    0.2696
%    0.4943
%   -1.4831
%   -1.0203
%   -0.4470
%    0.1097
%    1.1287
%   -0.2900
%    1.2616
%    0.4754
%    1.1741
%
% Cai Wingfield 2015-04
function element_list = elements(matrix_in)
    element_list = reshape(matrix_in, [], 1);
end%function
