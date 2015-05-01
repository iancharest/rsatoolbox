% itWasHeads = coinToss()
% itWasHeads = coinToss(probabilityOfHeads)
%
% Tosses a coin, and tells you whether it was heads or not. It's a 
% (pseudo-)fair coin unless otherwise specified.
%
% Good for using as conditionals.
%
%     % Randomly flip sign of x
%     if coinToss
%         x = (-1) * x;
%     end
% 
% CW 2015-04
function itWasHeads = coinToss(probabilityOfHeads)

    % If we don't specify otherwise, make it a fair coin.
    if nargin == 0
        probabilityOfHeads = 0.5;
    end

    % Put these edge cases in just in case the random roll comes up as
    % exactly 0 or 1, and we don't know now to treat < vs <=.
    if probabilityOfHeads == 0
        itWasHeads = false;
    elseif probabilityOfHeads == 1
        itWasHeads = true;
    else
        % rand is uniform over [0,1], so this should allow us to translate
        % from probability to threshold.
        itWasHeads = (rand <= probabilityOfHeads);
    end
end%function
