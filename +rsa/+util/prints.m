% function [stamped_message =] prints(FORMAT, A, ...)
%
% Sends the specified message to the output window, preceeded by a
% timestamp, and followed by a newline.  Optionally returns the string,
% minus the newline.
%
% Treats its arguments just as sprintf does, so formatting is possible.
%
% EXAMPLE USAGE
%
%     for i = 1:5
%         prints('Loop iteration number %d.', i);
%     end
%
% Produces something like:
%
%     [2015-03-24 17:29:40.968] Loop iteration number 1.
%     [2015-03-24 17:29:40.981] Loop iteration number 2.
%     [2015-03-24 17:29:40.994] Loop iteration number 3.
%     [2015-03-24 17:29:41.003] Loop iteration number 4.
%     [2015-03-24 17:29:41.011] Loop iteration number 5.
%
% See also WARNS
%
% Cai Wingfield 2015-03
function stamped_message = prints(varargin)

    % Get the current time asap.
    datestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
    
    % Apply the formatting as supplied in the argumets.
    message = sprintf(varargin{:});
    
    % Build the stamped string.
    stamped_message = ['[', datestamp, '] ', message];
    
    % Send it to the output window, including a trailing newline.
    fprintf([stamped_message, '\n']);
    
end%function
