% script adapted from toolbox for RFX on glm results

matlabpool open;

% all subjects 
subjectNames = { ...
    'meg08_0320', ...
    'meg08_0323', ...
    'meg08_0324', ...
    'meg08_0327', ...
    'meg08_0348', ...
    'meg08_0350', ...
    'meg08_0366', ...
    'meg08_0371', ...
    'meg08_0372', ...
    'meg08_0377', ...
    'meg08_0397', ...
    'meg08_0400', ...
    'meg08_0401', ...
    'meg08_0402' ...
    };

% loading predictors
predictors = { ...
    'Models-GTF1to30ms'; ...
    'VoiceModel'; ...
    'placeModel';
    'mannerModel';
    };

rootPath = '/imaging/iz01/NeuroLex/glm/';
significanceTestPermutations = 1000;

% for loop : predictor goes here
for predictor = 1:length(predictors);
    
    nSubjects = length(subjectNames);
    modelName = spacesToUnderscores(predictors{predictor});
    
    disp(modelName);
    
    fprintf('Permuting (RFX) ...');
    
    for subjectNumber = 1:nSubjects
        subject = subjectNames{subjectNumber};
        inputFilename = fullfile(rootPath, 'betaMaps', [subject '_' modelName]);
        
        MEGDataStcL = mne_read_stc_file1([inputFilename, '-lh.stc']);
        MEGDataStcR = mne_read_stc_file1([inputFilename, '-rh.stc']);
        MEGDataVolL = single(MEGDataStcL.data);
        MEGDataVolR = single(MEGDataStcR.data);
        all_rho(subjectNumber,:,:) = [MEGDataVolL; MEGDataVolR];
    end
    
    numberOfPermutation = significanceTestPermutations;
    [numberOfVertex, numberOfTimePoints] = size(MEGDataStcL.data);
    
    [h,p,ci,stats] = ttest(all_rho);
    t_value = squeeze(stats.tstat);
    t_lh = t_value(1:numberOfVertex,:);
    t_rh = t_value(numberOfVertex+1:numberOfVertex * 2,:);
    lh_Vol = MEGDataStcL;
    rh_Vol = MEGDataStcR;
    
    t_lh(isnan(t_lh)) = 0;
    t_rh(isnan(t_rh)) = 0;
    lh_Vol.data = t_lh;
    rh_Vol.data = t_rh;
    
    outputFilename = fullfile(rootPath, 'betaMaps', modelName, [modelName '_tMesh_allSubjects']);
    
    mne_write_stc_file1([outputFilename, '-lh.stc'], lh_Vol);
    mne_write_stc_file1([outputFilename, '-rh.stc'], rh_Vol);
    observed_Vol_L = lh_Vol;
    observed_Vol_R = rh_Vol;
    
    max_t_value = zeros(1,numberOfPermutation);
    
    parfor perm = 1:numberOfPermutation
        rho = all_rho;
        simulated_Vol_L = lh_Vol;
        simulated_Vol_R = rh_Vol;
        
        if mod(perm, floor(numberOfPermutation/40)) == 0, fprintf('\b.'); end%if
        
        for subjectNumber = 1:nSubjects
            %toss a coin
            toss = (rand(1) > 0.5)*2-1;
            rho(subjectNumber,:,:) = squeeze(all_rho(subjectNumber,:,:)).*toss;
        end
        [h,p,ci,stats] = ttest(rho);
        t_value = squeeze(stats.tstat);
        
        t_lh = t_value(1:numberOfVertex,:);
        t_rh = t_value(numberOfVertex+1:numberOfVertex * 2,:);
        simulated_Vol_L.data = t_lh;
        simulated_Vol_R.data = t_rh;
        
        
        outputFilename = fullfile(rootPath, 'betaMaps', modelName, ['perm-' num2str(perm) '_' modelName '_t_map']);
        mne_write_stc_file1([outputFilename, '-lh.stc'], simulated_Vol_L);
        mne_write_stc_file1([outputFilename, '-rh.stc'], simulated_Vol_R);
        
        max_t_value(perm) = max(max(t_value));
    end
    
    percent = 0.05; % Update from userOptions.primaryThreshold; IZ 03,12
    t_distribution = sort(max_t_value);
    
    vertex_level_threshold = t_distribution(ceil(size(t_distribution,2)*(1-percent)));
    
    fprintf('Writng results corrected at whole brain level using permutation but without using clustering method...\n');
    
    gotoDir(rootPath, 'Results');
    outputFileName_sig = fullfile(rootPath, 'Results', [modelName '_significant_vertex']);
    
    observed_Vol_L.data(observed_Vol_L.data<vertex_level_threshold) = 0;
    observed_Vol_R.data(observed_Vol_R.data<vertex_level_threshold) = 0;
    
    mne_write_stc_file1([outputFileName_sig, '-lh.stc'], observed_Vol_L);
    mne_write_stc_file1([outputFileName_sig, '-rh.stc'], observed_Vol_R);
    
    disp('Done!');

end

matlabpool close;