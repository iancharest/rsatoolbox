% [dp_min, dp_max] = toDataLimits(data_time_min, window_time_min, window_time_max, samplingRate)
%
% Converts the upper- and lower-limits of a time window to the upper- and
% lower-limits of data point indices.
%
%     data_time_min
%         The time in ms of the first data point.
%
%     window_time_min, window_time_max
%         The lower and upper limit, in ms, of the time window of interest.
%
%     samplingRate
%         The number of samples per second.
%
%     dp_min, dp_max
%         The lower and upper limit of the time window, in datapoint
%         indices.
%
% EXAMPLE USAGE
%
%     [~, allMEGData] = fiff_read_evoked(readPath);
%     samplingRate = allMEGData.info.sfreq/1000; % ms
%     tmin = double(allMEGData.evoked.first);
%     
%     [dp_min, dp_max] = toDataLimits(tmin, -200, 300, samplingRate)
%
% Based on code by IZ, 2012-09
function [dp_min, dp_max] = toDataLimits(data_time_min, window_time_min, window_time_max, samplingRate)
    dp_min = ((window_time_min - data_time_min)/samplingRate) + 1;
    dp_max = ((window_time_max - data_time_min)/samplingRate) + 1;
end%function
