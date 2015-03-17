% [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions)
%
% Given STCMetadata providing a picture of the data to be searched, and
% given the specification of the searchlight in miliseconds as set in
% userOptions, this function returns a struct which describes the
% specifications of the searchlight in datapoints.  Also returned will be
% slSTCMetadata, which is a new STC metadata struct, formatted for the
% searchlight specifications.
%
% STCMetadata must contain at least:
%     STCMetadata.tstep
%         The timestep of the data in seconds.
%     STCMetadata.tmin
%         The time of the first datapoint in seconds.
%     STCMetadata.tmax
%         The time of the last datapoint in seconds.
%
% slSTCMetadata will contain
%     slSTCMetadata.tstep
%         The timestep of the searchlight in seconds.
%     slSTCMetadata.tmin
%         The time of the first searchlight datapoint in seconds.
%     slSTCMetadata.tmax
%         The time of the last searchlight datapoint in seconds.
% ... and any other fields which STCMetadata had will also be copied.
%
% slSpec will contain:
%    slSpec.width
%        The width of the searchlight in datapoints.
%    slSpec.step
%        The step of the searchlight in datapoints.
%    slSpec.limits
%        The limits of the searchlight window of interest in datapoints.
%    slSpec.windowWidth
%        The number of datapoints in the searchlight window of interest.
%    slSpec.windowPositions
%        A (nPositions x 2)-matrix, where each row is the left and right
%        index of a position of the sliding window, in order.
%
% Based on code by IZ 2012
% Cai Wingfield 2015-03
function [slSpec, slSTCMetadata] = getSearchlightSpec(STCMetadata, userOptions)

    % TODO: this should also accept & return sensor-space data in place of
    % TODO: sensor-space data
    
    %% Common values

    % The timestep of the data in ms
    dataTimestep_data_ms = STCMetadata.tstep * 1000;
    
    % The time index of the first datapoint in ms
    firstPoint_data_ms = STCMetadata.tmin * 1000;
    
    % The number of timepoints in the data
    nTimepoints_data = (STCMetadata.tmax - STCMetadata.tmin) / STCMetadata.tstep;
    
    %% slSpec

    slSpec = struct();
    
    % The width in timepoints is...
    slSpec.width = ...
        ...% the width in ms...
        userOptions.temporalSearchlightWidth ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_data_ms;
    
    % The step in timepoints is...
    slSpec.step = ...
        ...% the timestep in ms...
        userOptions.temporalSearchlightTimestep ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_data_ms;
    
    slSpec.limits = [NaN, NaN];
    
    % The lower bound of the searchlight window of interest in timepoints
    % is ...
    slSpec.limits(1) = ...
        ...% the distance of the lower bound in ms from the start of the data...
        userOptions.temporalSearchlightLimits(1) - firstPoint_data_ms ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_data_ms;
    
    % The upper bound of the searchlight window of interest in timepoints
    % is...
    slSpec.limits(2) = ...
        ...% the distance of the upper bound in ms from the start of the data...
        userOptions.temporalSearchlightLimits(2) - firstPoint_data_ms ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_data_ms;
    
    % The width of the searchlight window is...
    slSpec.windowWidth = ...
        ...% the distance between the endpoints, plus one (fencepost).
        slSpec.limits(2) - slSpec.limits(1) + 1;
    
    % Calculate window positions
    
    % First window positions
    nextWindow = [slSpec.limits(1), slSpec.limits(1) + slSpec.width - 1];
    % List of window positions
    slSpec.windowPositions = nextWindow;
    % While we don't exceed the upper bound of the window...
    while nextWindow(2) + slSpec.step <= slSpec.limits(2)
        % ...move the window...
        nextWindow = nextWindow + slSpec.step;
        % ...and add it to the list of window positions.
        slSpec.windowPositions = [slSpec.windowPositions; nextWindow];
    end%while
    
    % Sanity check
    if slSpec.limits(1) < 1 || slSpec.limits(2) > nTimepoints_data
        error('getSearchlightSpec:InvalidSpec', 'Can''t produce a valid searchlight specification for this data.');
    end
    
    %% slSTCMetadata
    
    % Copy all existing fields from the input metadata struct
    slSTCMetadata = STCMetadata;
    
    % Need to convert from miliseconds as defined in the userOptions struct
    % to the seconds required by the STC metadata format.
    slSTCMetadata.tstep = userOptions.temporalSearchlightTimestep  / 1000;
    slSTCMetadata.tmin  = userOptions.temporalSearchlightLimits(1) / 1000;
    slSTCMetadata.tmax  = userOptions.temporalSearchlightLimits(2) / 1000;
        
end%function
