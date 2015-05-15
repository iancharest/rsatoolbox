% [e, s, worked] = pop(s)
%
% Given a list (stack) s, pop(s) returns the element at the top of the
% stack and the remaining stack.  The flag worked indicates whether the pop
% was successful.  For example, calling pop on an empty list will not be
% succesful, but if the list is not empty, it should be.
%
% CW 2015-05
function [e, s, worked] = pop(s)
    % The stack will be a column vector.
    s = s(:);
    if isempty(s)
        e = nan;
        worked = false;
        return;
    elseif numel(s) == 1
        e = s(1);
        s = [];
        worked = true;
        return;
    else
       e = s(1);
       s = s(2:end);
       worked = true;
       return;
    end
end%function
