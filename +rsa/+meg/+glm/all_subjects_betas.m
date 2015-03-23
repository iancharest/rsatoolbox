% all subjects individual model bets calculation

% loading responses
path = '/imaging/iz01/NeuroLex/cochleagram_HD/RDMs/';

%%%%%%%%%%%%%%%%
% all subjects %
%%%%%%%%%%%%%%%%
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


output = mne_read_stc_file1(['/imaging/iz01/NeuroLex/cochleagram_HD/Results/lexPro_tonotopy_1to30ms_HD_Gamma30ms_significant_vertex-lh.stc']);
output.tmin = -0.05;
output.tstep = 0.01;

% scale for coloring
% scale = [15; 5; 5; -15]; % cochlear, voice, place, manner
% models = {'cochlear', 'voice', 'place', 'manner'};


for i = 1:length(subjectNames) % all subjects
    
    disp(subjectNames{i});
 
    for j = 1:2 % both hemispheres
        
        if j==1 chi = 'L'; disp( '- left hemisphere');
        else chi = 'R'; disp(' - right hemisphere'); end
        
        betas = zeros(10242,22);
        load([path 'searchlightRDMs_masked_' subjectNames{i} '-' lower(chi) 'h.mat']);
        
        % applying glm
        vertices = fieldnames(searchlightRDMs);
        
        for predictor = 1:length(predictors) % all models
            disp([ ' -- ' predictors{predictor}]);
            load(predictors{predictor});
            eval(['model =' strtok(predictors{predictor},'-') '.RDM;'] );
            
            for v =1:length(vertices)
                thisVertex = vertices{v};
                [temp vertex_num] = strtok(thisVertex,'_'); vertex_num = str2num(strtok(vertex_num,'_'));
                nTimePoints = fieldnames(searchlightRDMs.(thisVertex));
                
                for t = 1:length(nTimePoints)
                    thisTime = nTimePoints{t};
                    
                    [b,dev,stats] = glmfit(vectorizeRDM(model)',...
                        vectorizeRDM(searchlightRDMs.(thisVertex).(thisTime).RDM)','normal');
                    
                    betas(vertex_num,t) = b(2,1);
                    
%                     disp([thisVertex ' ' thisTime ' ' num2str(b{v,t}(2,1),4)]);
                    
                end % time

            end % vertex
            
            % output settings
            
            output.data = betas;
            
            mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/betaMaps/' subjectNames{i} '_' predictors{predictor} '-' lower(chi) 'h.stc'], output);
        end % predictor
        

    end % hemisphere
    disp(['Done for subject: ' subjectNames{i}]);

end % subject
    
% save('glmResults','output','b','scale','models','dev','stats');
 
 