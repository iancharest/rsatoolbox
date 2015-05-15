% Given an element e and an array (stack) s, 
% s = push(e, s) pushes e into the first position of s.
%
% CW 2015-05
function s = push(e, s)
    % The stack will be a column vector.
    s = [e; s(:)];
end%function
