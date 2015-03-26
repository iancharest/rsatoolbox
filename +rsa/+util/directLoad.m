% [x, y, ...] = directLoad(pathToLoad, {'x', 'y', ...})
%
% Whereas the built-in load() function will only load variables into a
% struct with variable names for fields, directLoad will load the
% variables' contents directly.
%
% So
%     x = load(path, 'x');
%     x = x.x;
% is equivalent to
%     x = directLoad(path, 'x');
%
% And
%     a = load(path, {'x', 'y'});
%     x = a.x;
%     y = a.y;
% is equivalent to
%     [x, y] = directLoad(path, {'x', 'y'});
%
% Cai Wingfield 2015-03
function [varargout] = directLoad(pathToLoad, varNames)
    % Make sure we have a cell array of variable names to load
    if ~iscell(varNames)
        varNames = {varNames};
    end
    
    nVars = numel(varNames);
    loadedVars = load(pathToLoad, varNames{:});
    varargout = cell(1, nVars);
    for var_i = 1:nVars
        varName = varNames{var_i};
        varargout{var_i} = loadedVars.(varName);
    end
end%function
