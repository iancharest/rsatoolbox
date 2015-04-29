% [x, y, ...] = directLoad(pathToLoad[, {'x', 'y', ...}])
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
% Alternatively, if there is only one variable saved in the file, its name
% can be ommitted.
%
%     a = directLoad('path/to/file/containing/one-variable.mat');
%
% However, if this form is used when there is more than one variable saved
% in the file, an error will be thrown.
%
% Cai Wingfield 2015-03--2015-04
function [varargout] = directLoad(pathToLoad, varNames)
    
    % Make sure we have a cell array of variable names to load
    if nargin == 1
        varNames = {};
    elseif ~iscell(varNames)
        varNames = {varNames};
    end
    
    % If there are no specified varnames, we load all variables in the
    % file.  We require that there is only one such variable if this is the
    % case.
    if isempty(varNames)
        loaded_variables = load(pathToLoad);
        loaded_variable_names = fieldnames(loaded_variables);
        
        if numel(loaded_variable_names) ~= 1
           error('The file "%s" contains multiple variable names. Please specify which you want to load.', pathToLoad); 
        else
            varargout{1} = loaded_variables.(loaded_variable_names{1});
        end%if
        
    % If there are specified variable names to load, we load them.
    else
        nVars = numel(varNames);
        loaded_variables = load(pathToLoad, varNames{:});
        varargout = cell(1, nVars);
        for var_i = 1:nVars
            varName = varNames{var_i};
            varargout{var_i} = loaded_variables.(varName);
        end%for
    end%if
end%function
