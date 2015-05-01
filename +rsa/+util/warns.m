% function [stamped_message =] warns(FORMAT, A, ...)
%
% Sends the specified message as a warning to the output window, preceeded 
% by a timestamp, and followed by a newline.  Optionally returns the
% string, minus the newline.
%
% Treats its arguments just as sprintf does, so formatting is possible.
%
% Can be suppressed with `warning off`.
%
% EXAMPLE USAGE
%
%     for i = 1:5
%         warns('Loop iteration number %d.', i);
%     end
%
% Produces something like:
%
%     [2015-03-24 17:29:40.968] Warning: Loop iteration number 1.
%     [2015-03-24 17:29:40.981] Warning: Loop iteration number 2.
%     [2015-03-24 17:29:40.994] Warning: Loop iteration number 3.
%     [2015-03-24 17:29:41.003] Warning: Loop iteration number 4.
%     [2015-03-24 17:29:41.011] Warning: Loop iteration number 5.
%
% See also rsa.util.PRINTS rsa.util.ERRORS
%
% Cai Wingfield 2015-03
function stamped_message = warns(varargin)

    % Get the current time asap.
    datestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
    
    % Apply the formatting as supplied in the argumets.
    message = sprintf(varargin{:});
    
    % Build the stamped string.
    stamped_message = ['[', datestamp, '] Warning: ', message];
    
    % If warnings are all turned off, we suppress sending to stdout.
    warning_state = warning;
    if any(ismember({warning_state.state}, 'on'))
        % Send it to the output window, including a trailing newline.
        fprintf(['[\b' stamped_message, ']\b\n']);
    end
    
end%function
