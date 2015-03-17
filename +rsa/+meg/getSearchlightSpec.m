% slSpec = getSearchlightSpec(STCMetadata, userOptions)
%
% Given STCMetadata providing a picture of the data to be searched, and
% given the specification of the searchlight in miliseconds as set in
% userOptions, this function returns a struct which describes the
% specifications of the searchlight in datapoints.
%
% STCMetadata must contain at least:
%     STCMetadata.tstep
%         The timestep of the data in seconds.
%     STCMetadata.tmin
%         The time of the first datapoint in seconds.
%     STCMetadata.tmax
%         The time of the last datapoint in seconds.
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
function slSpec = getSearchlightSpec(STCMetadata, userOptions)

    % TODO: this should also accept sensor-space data in place of
    % TODO: sensor-space data

    % The timestep of the data in ms
    dataTimestep_ms = STCMetadata.tstep * 1000;
    
    % The time index of the first datapoint in ms
    firstPoint_ms = STCMetadata.tmin * 1000;
    
    % The number of timepoints in the data
    nTimepoints = (STCMetadata.tmax - STCMetadata.tmin) / STCMetadata.tstep;

    slSpec = struct();
    
    % The width in timepoints is...
    slSpec.width = ...
        ...% the width in ms...
        userOptions.temporalSearchlightWidth ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_ms;
    
    % The step in timepoints is...
    slSpec.step = ...
        ...% the timestep in ms...
        userOptions.temporalSearchlightTimestep ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_ms;
    
    slSpec.limits = [NaN, NaN];
    
    % The lower bound of the searchlight window of interest in timepoints
    % is ...
    slSpec.limits(1) = ...
        ...% the distance of the lower bound in ms from the start of the data...
        userOptions.temporalSearchlightLimits(1) - firstPoint_ms ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_ms;
    
    % The upper bound of the searchlight window of interest in timepoints
    % is...
    slSpec.limits(2) = ...
        ...% the distance of the upper bound in ms from the start of the data...
        userOptions.temporalSearchlightLimits(2) - firstPoint_ms ...
        ...% divided by the timestep of the data in ms...
        / dataTimestep_ms;
    
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
    if slSpec.limits(1) < 1 || slSpec.limits(2) > nTimepoints
        error('getSearchlightSpec:InvalidSpec', 'Can''t produce a valid searchlight specification for this data.');
    end
        
end%function
