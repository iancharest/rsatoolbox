% time_to_give_feedback = feedback_throttle(target_feedback_count, current_loop_iteration, total_loop_iterations)
%
% Allows throtting of feedback in loops.
%
% Use it like this:
%
%     loop_count = 110;
%     n_feedbacks = 4;
%     for i = 1:loop_count
%         % ...
%         % Do something that takes ages
%         % ..
%         if feedback_throttle(n_feedbacks, i, loop_count)
%             rsa.util.prints( ...
%                 'Now %2.0f%% complete.', rsa.stat.percent(i, loop_count);
%         end
%     end
%
% Which results in:
%
%     [2015-05-01 14:36:02.591] Now 24% complete.
%     [2015-05-01 14:49:25.833] Now 49% complete.
%     [2015-05-01 15:03:11.112] Now 73% complete.
%     [2015-05-01 14:12:57.320] Now 98% complete.
%
% See also rsa.stat.percent rsa.util.prints
%
% CW 2015-05
function time_to_give_feedback = feedback_throttle(target_feedback_count, current_loop_iteration, total_loop_iterations)
    increments_per_feedback = floor(total_loop_iterations/target_feedback_count);
    time_to_give_feedback = (mod(current_loop_iteration, increments_per_feedback) == 0);
end%function
