function averageRDMPaths = averageSearchlightRDMs(RDMPaths, userOptions)

    import rsa.util.*
    
    % Paths
    file_i = 1;
    for chi = 'LR'
        averageRDMPaths.(chi) = fullfile(userOptions.rootPath, 'RDMs', ['average_', lower(chi), 'h.mat']);
        promptOptions.checkFiles(file_i).address = averageRDMPaths.(chi);
        file_i = file_i + 1;
    end
    
    promptOptions.functionCaller = 'MEGDataPreparation_source';
    promptOptions.defaultResponse = 'S';
    
    overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
    if overwriteFlag

        nSubjects = numel(userOptions.subjectNames);
        for chi = 'LR'

            nVertices = NaN;
            nTimepoints = NaN;

            for subject_i = 1:nSubjects

                this_subject_name = userOptions.subjectNames{subject_i};

                prints('Loading searchlight RDMs for subject %s (%d/%d) %sh...', this_subject_name, subject_i, nSubjects, lower(chi));
                this_subject_slRDMs = directLoad(RDMPaths(subject_i).(chi), 'searchlightRDMs');

                % For the first subject, we initialise the average and the
                % nan-counter with some sizes.
                if subject_i == 1
                    [nVertices, nTimepoints] = size(this_subject_slRDMs);
                    average_slRDMs(1:nVertices, 1:nTimepoints) = struct('RDM', zeros(size(this_subject_slRDMs(1,1).RDM)));
                    nan_counts(1:nVertices, 1:nTimepoints) = struct('mask', zeros(size(this_subject_slRDMs(1,1).RDM)));
                end

                prints('Adding RDMs at all vertices and timepoints...');

                parfor t = 1:nTimepoints
                    for v = 1:nVertices
                        nan_locations = isnan(this_subject_slRDMs(v, t).RDM);
                        this_subject_slRDMs(v, t).RDM(nan_locations) = 0;
                        average_slRDMs(v, t).RDM = average_slRDMs(v, t).RDM + this_subject_slRDMs(v, t).RDM;
                        nan_counts(v, t).mask = nan_counts(v, t).mask + nan_locations;
                    end%for:t
                end%for:v
            end%for:subject

            prints('Averaging RDMs at all vertices...');

            % replace nan counts by non-nan counts
            parfor t = 1:nTimepoints
                for v = 1:nVertices
                    non_nan_counts = nSubjects - nan_counts(v, t).mask;
                    average_slRDMs(v, t).RDM = average_slRDMs(v, t).RDM ./ non_nan_counts;
                end
            end

            prints('Saving average searchlight RDMs to "%s"...', averageRDMPaths.(chi));
            save('-v7.3', averageRDMPaths.(chi), 'average_slRDMs');

        end%for:chi

    else
        prints('Average RDMs already calculated.  Skipping...');
    end
end%function
