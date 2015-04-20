% s = sumall(M)
%
% Sums all elements in an n-dimensional array.
%
% EG:
% >> a = randn(2,2,2,2);
% >> sumall(a)
%  ans =
%      0.9533
%
% Cai Wingfield 2015-04
function s = sumall(M)
    import rsa.util.*

    s = sum(elements(M));
end%function
