% function errors(FORMAT, A, ...)
%
% Sends the specified message to the output window as an error, preceeded
% by a timestamp, and followed by a newline.
%
% Treats its arguments just as sprintf does, so formatting is possible.
%
% EXAMPLE USAGE
%
%     problem = 'there was a problem';
%     errors('Terminated with the message "%s".', problem);
%
% Produces something like:
%
%     Error using rsa.util.errors (line 36)
%     [2015-05-01 11:27:30.769] Terminated with the message "there was a problem".
%
% See also rsa.util.WARNS rsa.util.PRINTS
%
% Cai Wingfield 2015-04
function errors(varargin)

    % Get the current time asap.
    datestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
    
    % Apply the formatting as supplied in the argumets.
    message = sprintf(varargin{:});
    
    % Build the stamped string.
    stamped_message = ['[', datestamp, '] ', message];
    
    % Send it to the output window, including a trailing newline.
    error(stamped_message);
    
end%function
