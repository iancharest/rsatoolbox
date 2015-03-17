% This function uses userOptions to first generate masked source meshes for
% current time window and calculated group statistics across sliding time
% windows. Generates plots and saves all data RDMs and computed statistics
% as excel spreadsheet.
%
% Script can perform two kinds of group statistics:
%   Random Effects Analysis: This will generate an uncorrected t-map and a
%       corrected t-map thresholded based on cluster statistics. Each cluster
%       will have a mass and a corresponding p-value.
%
%   Fixed Effects Analysis: This will generate correlation values with
%       corresponding p values.

% Based on scripts by Li Su
% Written by Isma Zulfiqar 12-12
% Updated by IZ 03/13 Now "only" computes data RDMs using sliding time window
% updated by FJ 05/12/13 fixing cases senitive issues
% Updated by FJ 03/14

function ROI_slidingTimeWindow(userOptions, model)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
import rsa.par.*

returnHere = pwd; % We'll come back here later
modelName = model.name;
if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

readPathL = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', userOptions.betaCorrespondence(1, 1).identifier, '[[subjectName]]', userOptions.subjectNames{1}, '[[LR]]', 'l');
MEGDataStcL = mne_read_stc_file1(readPathL);

output_path = fullfile(userOptions.rootPath, 'RDMs');
output_path_Meshes= fullfile(userOptions.rootPath, 'RDMs','Mesh');
mkdir(output_path,'Mesh');
RDMs_filename = [userOptions.analysisName '_' modelName '_dataRDMs_sliding_time_window.mat'];
promptOptions.functionCaller = 'ROI_slidingTimeWindow';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(output_path, RDMs_filename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag

    nSubjects = numel(userOptions.subjectNames);
    
    sourceMesh_load = struct('L',[], 'R',[]);
 
    parfor subject =1:nSubjects
       disp(['Reading source meshes for ' num2str(subject) '... ']);
       sourceMesh_load = MEGDataPreparation_source (subject,betaCorrespondence, userOptions);
       parsave(fullfile(output_path_Meshes, userOptions.subjectNames{subject}), sourceMesh_load);
    end
    
    prints('Done with Reading Source Meshes ');
        
    for subject =1:nSubjects

        % Data preparation %
        thisSubject = userOptions.subjectNames{subject};
        % RDM Calculation %
        fprintf('Computing data RDMs...');
        nMasks = numel(userOptions.maskNames);
        localOptions = userOptions;
        localOptions.subjectNames = {thisSubject};
        localOptions.nSubjects = 1;
        
        sourceMesh_load = load (fullfile(output_path_Meshes, userOptions.subjectNames{subject}));
        
        % TODO: make sure this uses getSearchlightSpec properly
        
        for mask=1:nMasks
            thisMask = userOptions.maskNames{mask};
            disp(thisMask);
            currentTimeWindow = userOptions.maskTimeWindows{mask};

            timeWindow=1;
            if currentTimeWindow(1) < MEGDataStcL.tmin*1000, currentTimeWindow(1) = MEGDataStcL.tmin*1000; end
            for time = currentTimeWindow(1):userOptions.temporalSearchlightTimestep:currentTimeWindow(2)-userOptions.temporalSearchlightWidth
                disp([num2str(time) ' ms:'])
                
                localOptions.maskNames = {thisMask};
                localOptions.maskTimeWindows = {[time time+userOptions.temporalSearchlightWidth]};
                tw{timeWindow} = num2str(localOptions.maskTimeWindows{1});
                
                lowerTimeLimit = time;
                upperTimeLimit = time + userOptions.temporalSearchlightWidth; % i hope this makes sense!
                
                % getting the starting data point by comparing starting points of searchlight and input data
                differenceInms = lowerTimeLimit - MEGDataStcL.tmin*1000;
                startingDataPoint = 1 + floor((differenceInms / (localOptions.STCmetaData.tstep*1000)));
                
                % getting last data point for searchlight limits
                differenceInms = norm(upperTimeLimit - lowerTimeLimit);
                lastDataPoint = startingDataPoint + floor((differenceInms / (MEGDataStcL.tstep*1000))) -1;
                
                if lastDataPoint <= userOptions.maskTimetoDataPoints.(dashToUnderscores(thisMask))(2)+ceil(userOptions.temporalSearchlightWidth/userOptions.temporalSearchlightTimestep)
                    localOptions.maskTimetoDataPoints.(dashToUnderscores(thisMask)) = [startingDataPoint lastDataPoint];
                    
                    indexMasks = MEGMaskPreparation_source(localOptions);
                    
                    maskedMeshes = MEGDataMasking_source(sourceMesh_load, indexMasks, localOptions.betaCorrespondence, localOptions);
                    
                    %% RDM calculation %%
                    
                    RDMs = constructRDMs(maskedMeshes, localOptions.betaCorrespondence, localOptions);
                    % RDMs = averageRDMs_subjectSession(RDMs, 'session');
                    
                    RDMs.RDM = vectorizeRDM(RDMs.RDM);
                    RDMs.name = [thisMask ' | ' num2str(time) 'ms | ' thisSubject];
                    allRDMs(mask,timeWindow,subject,1:localOptions.nSessions) = RDMs;
                    
                    timeWindow = timeWindow+1;                    
                else
                    disp('Time point does not exist in the data.');
                    break;
                end
            end % for time
            
        end % for:mask
        
    end
    
    fprintf('Saving all data RDMs... ');
    save ('-v7.3',fullfile(output_path, RDMs_filename), 'allRDMs');
    disp('Done!');
    rmdir(output_path_Meshes, 's');
else
    fprintf('Data RDMs have already been computed, skip....\n');
    
end

