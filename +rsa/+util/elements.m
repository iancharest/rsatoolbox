% ms = elements(M)
%
% Matlab won't let you use (:) inline as in
%     sum(x)(:)
% so this functon lets that happen:
%     elements(sum(x))
%
% CW 2015-08
function ms = elements(M)
    ms = M(:);
end%function
