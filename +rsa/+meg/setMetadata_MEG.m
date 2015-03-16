% userOptions = setMetadata_MEG(Models, userOptions)
%
% Models - an array of structs
% userOptions - struct of user preferences, modified and returned.
%
% This function sets any missing parameters in userOptions to the default
% values. If no default values can be set, ask the user to set it and show
% an error message. It also reads in the MEG data to find out some meta
% data related to the experiment.
%
% Li Su 3-2012
% Update: Isma Zulfiqar 9-2012 fixed time issues for searchlight all brain
% Update: IZ 11-2012 time support for ROI maks

function userOptions = setMetadata_MEG(Models, userOptions)

% loading one subject's STC/FIFFfile in order to find out its meta data.
tempBetas = betaCorrespondence();
userOptions.betaCorrespondence = tempBetas;

usingMasks = ~isempty(userOptions.maskNames);

if userOptions.sensorLevelAnalysis
    
    readPath = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', tempBetas(1,1).identifier, '[[subjectName]]', userOptions.subjectNames{1});
    
    [ignore, allMEGData] = fiff_read_evoked(readPath);
    
    samplingRate = allMEGData.info.sfreq/1000; % ms
    
    tmin = double(allMEGData.evoked.first); tmax = double(allMEGData.evoked.last);
    
    % TODO: userOptions.searchlight isn't in use any more. This adjustment
    % TODO: should be done once, depending on the recipe being used.
    if userOptions.searchlight
        time = userOptions.temporalSearchlightLimits;
        
        if time(1) < tmin
            disp(['Lower time limit is out of bounds. Using minimum value from fiff file instead. tmin = ' num2str(tmin) ' ms']);
            userOptions.temporalSearchlightLimits(1) = tmin;
        end
        if time(2) > tmax
            disp(['Upper time limit is out of bounds. Using maximum value from fiff file instead. tmax = ' num2str(tmax) ' ms']);
            userOptions.temporalSearchlightLimits(2) = tmax;
        end
        
        userOptions.toDataPoints = [1+((userOptions.temporalSearchlightLimits(1) - tmin)/samplingRate) 1+((userOptions.temporalSearchlightLimits(2)-tmin)/samplingRate)];
        
    else
        time = userOptions.maskSpec.timeWindow;
        
        if time(1) < tmin
            disp(['Lower time limit is out of bounds. Using minimum value from fiff file instead. tmin = ' num2str(tmin) ' ms']);
            userOptions.maskSpec.timeWindow(1) = tmin;
        end
        if time(2) > tmax
            disp(['Upper time limit is out of bounds. Using maximum value from fiff file instead. tmax = ' num2str(tmax) ' ms']);
            userOptions.maskSpec.timeWindow(2) = tmax;
        end
        
        userOptions.maskSpec.toDataPoints = [1+((userOptions.maskSpec.timeWindow(1) - tmin)/samplingRate) 1+((userOptions.maskSpec.timeWindow(2)-tmin)/samplingRate)];
        
    end
    
else % source level analysis
    
    % We can read just the left hemisphere, as we just want the metadata,
    % and we assume both hemispheres will be the same
    readPathL = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', tempBetas(1, 1).identifier, '[[subjectName]]', userOptions.subjectNames{1}, '[[LR]]', 'l');
    MEGDataStcL = mne_read_stc_file1(readPathL);
    MEGDataVolL = single(MEGDataStcL.data);
    
    % TODO: this is bad
    modelNumber = userOptions.modelNumber;
    
    % Todo: this is bad
    userOptions.nConditions = size(squareRDM(Models(modelNumber).RDM), 1);
    
    % ============= setting time parameters for output file ================= %
    userOptions.STCmetaData.tmin = MEGDataStcL.tmin; % - ...
    % (MEGDataStcL.tmin - (userOptions.temporalSearchlightLimits(1) /1000)); % in seconds
    userOptions.STCmetaData.vertices = 1:userOptions.targetResolution;
    userOptions.STCmetaData.tstep = MEGDataStcL.tstep;
    
    % time values are converted to their data point equivalents here.
    % this depends on the size of recordings input
    
    % ================= converting from mseconds to datapoints ============== %
    
    [vertices, totalDataPoints] = size(MEGDataStcL.data);
    totalTimeInMs = totalDataPoints * MEGDataStcL.tstep*1000; % calculated from time step and total data points
    timeAdjusted = totalTimeInMs - norm(MEGDataStcL.tmin * 1000, 2); % last point of data from tmin and total time
    
    %============================= input checks ============================= %
    % ====== comparing search light resolution to the time step of data ===== %
    
    % TODO: userOptions.searchlight isn't in use any more. This adjustment
    % TODO: should be done once, depending on the recipe being used.
    % time step
    if userOptions.searchlight % added 03/12 IZ
        if MEGDataStcL.tstep * 1000 * userOptions.temporalDownsampleRate > userOptions.temporalSearchlightResolution
            error('Error: The input time resolution of search light cannot be applied. The time resolution for data is lower.');
        else
            userOptions.STCmetaData.tstep = userOptions.temporalSearchlightResolution * ...
                userOptions.temporalDownsampleRate / 1000;  % in s
        end
    end
    
    % searchlight width in relation to downsample rate
    if userOptions.temporalSearchlightWidth < userOptions.STCmetaData.tstep * 1000
        disp('Warning: The searchlight width is less than data rate after downsampling');
        disp(['Setting temporal searchlight width equal to minimum timestep: ', num2str(userOptions.STCmetaData.tstep*1000), 'ms']);
        userOptions.temporalSearchlightWidth = userOptions.STCmetaData.tstep*1000;
    end
    
    %% MASKING TIME WINDOWS %%
    % TODO: userOptions.searchlight isn't in use any more. This adjustment
    % TODO: should be done once, depending on the recipe being used.
    if usingMasks && not(userOptions.searchlight) % sliding time window and roi analysis
        
        nMasks = numel(userOptions.maskNames);
        for mask_i = 1:nMasks
            
            % Which mask is this?
            thisMaskName = dashToUnderscores(userOptions.maskNames{mask_i});
            
            time = cell2mat(userOptions.maskTimeWindows(mask_i));
            lowerTimeLimit = time(1);
            upperTimeLimit = time(2);
            
            % getting the starting data point by comparing starting points of searchlight and input data
            differenceInms = lowerTimeLimit - MEGDataStcL.tmin*1000;
            startingDataPoint = 1 + floor((differenceInms / (MEGDataStcL.tstep*1000)));
            
            % getting last data point for searchlight limits
            differenceInms = norm(upperTimeLimit - lowerTimeLimit);
            lastDataPoint = startingDataPoint + floor((differenceInms / (MEGDataStcL.tstep*1000)));
            
            % parameters for output file
            userOptions.STCmetaData.tmin = lowerTimeLimit /1000; % in seconds
            
            % checks onmask timing information
            if lowerTimeLimit < MEGDataStcL.tmin*1000
                disp(['Warning: The searchlight for mask: ', thisMaskName, ' is attempting to access time point not in data.'] );
                disp(strcat(' >> Using minimum time point value instead... ', int2str(MEGDataStcL.tmin*1000), ' ms'));
                startingDataPoint = 1;
                userOptions.STCmetaData.tmin = MEGDataStcL.tmin;
            end
            
            % TODO: slidingTimeWindow isn't used anymore.  This adjustment
            % TODO: should be done in the relevant recipes.
            if userOptions.slidingTimeWindow % added by IZ 11-12
                if timeAdjusted-userOptions.temporalSearchlightWidth <= upperTimeLimit
                    disp(['Warning: The sliding window for mask: ', thisMaskName, ' over-runs the data sample. ']);
                    disp(strcat('>> Using maximum time point, with adjusted factor of source searchlight radius, from data instead... ', int2str(timeAdjusted-userOptions.sourceSearchlightRadius), ' ms'));
                    lastDataPoint = totalDataPoints - ...
                        ceil(userOptions.temporalSearchlightWidth/userOptions.temporalSearchlightResolution);
                end
            else
                if timeAdjusted <= upperTimeLimit
                    disp(['Warning: The searchlight for mask: ', thisMaskName, ' over-runs the data sample. ']);
                    disp(strcat('>> Using maximum time point from data instead... ', int2str(timeAdjusted-userOptions.sourceSearchlightRadius), ' ms'));
                    lastDataPoint = totalDataPoints;
                end
            end
            
            % TODO: Don't store this in the userOptions
            userOptions.maskTimetoDataPoints.(thisMaskName) = [startingDataPoint lastDataPoint];
        end
        
    % TODO: userOptions.searchlight isn't in use any more. This adjustment
    % TODO: should be done once, depending on the recipe being used.
    elseif userOptions.searchlight % all brain searchlight
        
        userOptions.STCmetaData.tmin = userOptions.temporalSearchlightLimits(1)/1000;
        
        % getting the starting data point by comparing starting points of searchlight and input data
        differenceInms = userOptions.temporalSearchlightLimits(1) - MEGDataStcL.tmin*1000;
        startingDataPoint = 1 + floor((differenceInms / (MEGDataStcL.tstep*1000)));
        
        % getting last data point for searchlight limits
        differenceInms = norm(userOptions.temporalSearchlightLimits(2) - userOptions.temporalSearchlightLimits(1));
        lastDataPoint = startingDataPoint + floor((differenceInms / (MEGDataStcL.tstep*1000)));
        
        % lower temporal searchlight limit
        if userOptions.temporalSearchlightLimits(1) < MEGDataStcL.tmin*1000
            disp('Warning: The searchlight is attempting to access time point not in data.' );
            disp(strcat(' >> Using minimum time point value instead... ', int2str(MEGDataStcL.tmin*1000), ' ms'));
            startingDataPoint = 1;
            userOptions.STCmetaData.tmin = MEGDataStcL.tmin;
        end
        
        % upper temporal searchlight limit
        if timeAdjusted-userOptions.temporalSearchlightWidth <= userOptions.temporalSearchlightLimits(2)
            disp('Warning: The search light over-runs the data sample. ');
            disp(strcat('>> Using maximum time point, with adjusted factor of source searchlight radius, from data instead... ', int2str(timeAdjusted-userOptions.sourceSearchlightRadius), ' ms'));
            userOptions.temporalSearchlightLimits(2) = timeAdjusted;
            lastDataPoint = totalDataPoints - ...
                (userOptions.temporalSearchlightWidth/userOptions.temporalSearchlightResolution);
        end
        
        userOptions.dataPointsSearchlightLimits = [startingDataPoint lastDataPoint];
        
    end
    
    % ========================== new fields ================================ %
    userOptions.toDataPoints = totalDataPoints/(userOptions.temporalDownsampleRate*totalTimeInMs);
    
end
