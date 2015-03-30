toolboxRoot = '/imaging/iz01/test/toolbox/devel/toolbox'; addpath(genpath(toolboxRoot));
% loading response
path = '/imaging/ls02/NeuroLex/cochleagram_HD/RDMs/';
chi = 'R';
nbins = 8;

nSubjects = 1;

% loading predictors
load('dynamic_models_8bins.mat');

load('frqs.mat');
frq_bins = frqs;

temppath = '/imaging/iz01/NeuroLex/cochleagram_HD/Results/';

output = mne_read_stc_file1([temppath 'lexPro_tonotopy_1to30ms_HD_Gamma30ms_significant_vertex-lh.stc']);
output.tmin = -0.05;
output.tstep = 0.01;

scale = [-25; -12.5; -6.2; -3.1; 3.1; 6.2; 12.5; 25;]; 
models = {'bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8'};

betas_glm_amp = zeros(10242,22,nSubjects);
betas_glm_scale = zeros(10242,22,nSubjects);

betas_ksdensity_amp = zeros(10242,22,nSubjects);
betas_ksdensity_cf = zeros(10242,22,nSubjects);

betas_normfit_sd = zeros(10242,22,nSubjects);
betas_normfit_cf = zeros(10242,22,nSubjects);

for i = 1:nSubjects
%     disp(subjectNames{i});
    disp(' -Loading data RDM...')
    load([path 'averaged_searchlightRDMs_masked_' lower(chi) 'h.mat']);
    
    % applying glm
    vertices = fieldnames(averageSubjectRDMs.(chi));

    for v =1:length(vertices)
        
        thisVertex = vertices{v};
        [temp vertex_num] = strtok(thisVertex,'_'); vertex_num = str2num(strtok(vertex_num,'_'));
        nTimePoints = fieldnames(averageSubjectRDMs.(chi).(thisVertex));
        
        disp([' --vertex: ' num2str(vertex_num)]);
        count=1;
        for t = 1:length(nTimePoints)
            thisTime = nTimePoints{t};
            disp([' ---time point: ' thisTime]);
            
            all_models =[];
            for m=1:nbins
                all_models = [all_models ; vectorizeRDM(allmodels{m,t}.RDM);];
            end
            
            [b,dev,stats] = glmfit(all_models',...
                vectorizeRDM(averageSubjectRDMs.(chi).(thisVertex).(thisTime).RDM),'normal');
            
            
            which_model = find(b==max(b(2:end,:))) - 1; % 1's are added by glm as predictor
            
            betas_glm_amp(vertex_num, t, i) = max(b(2:end,:));
            betas_glm_scale(vertex_num,t,i) = scale(which_model);
            
            disp([' ---maximum beta: ' num2str(max(b(2:end,:)))]);
            disp([' ---model with maximum beta: ' models{which_model}]);
                        
            base = ceil(b(2:end,:) * 10000);
            base(base < 0) = 0;
            
            % curve fitting
            h=[];
            h(1:base(1)) = 1; ind =[1];
            for index = 2:length(base)   
                h(sum(base(ind)) + 1:sum(base([ind; index]))) = index;
                ind = [ind; index];
            end
            if ~isempty(h)
                [f,x] = ksdensity(h);
%                 subplot(4,1,count)
                plot(x,f*300,'color','r');
                hold on
                hist(h,length(models));
                
                bin = floor(x(f==max(f))); % bin
                bin = bin(1);
                if bin==0
                    bin=1;
                    cf=1;
                else
                    cf = floor((x(f==max(f)) - bin) * nbins); % cf
                    cf = cf(1);
                    if cf==0, cf=1; end
                end
                betas_ksdensity_amp(vertex_num,t,i) = max(f)*1000;
                betas_ksdensity_cf(vertex_num,t,i) = frq_bins(cf,bin);
                
                disp([' ---ksdensity cf:' num2str(cf) ' bin:' num2str(bin)]);
                
                % normfit
                [mu, sig] = normfit(h);
                bin = floor(mu);
                cf = floor((mu-bin) * nbins); if cf==0, cf=1; end
                betas_normfit_sd(vertex_num,t,i) = sig;
                betas_normfit_cf(vertex_num,t,i) = frq_bins(cf,bin);
                
                disp([' ---normfit cf:' num2str(cf) ' bin:' num2str(bin)]);
                
                count = count+1;
                if mod(t,4)==0,         title(['Vertex: ' num2str(vertex_num)]);
                    %                 saveas(figure(1),['/imaging/iz01/NeuroLex/glm/ksdensityFigures/' chi '_' num2str(subjectNames{i}) '-v' num2str(vertex_num) '-t' num2str(t)],'fig');
                    close all; figure; count = 1;
                end
            else
                betas_ksdensity_amp(vertex_num,t,i) = NaN;
                betas_ksdensity_cf(vertex_num,t,i) = NaN;
                
                betas_normfit_sd(vertex_num,t,i) = NaN;
                betas_normfit_cf(vertex_num,t,i) = NaN;
                
            end % if h not empty
        end % time
    end % vertex
    output.data = squeeze(betas_glm_amp(:,:,i));
    mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/dynamicModels/eightbins/avg_glm_beta-' lower(chi) 'h.stc'], output);
    
    output.data = squeeze(betas_glm_scale(:,:,i));
    mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/dynamicModels/eightbins/avg_glm_scale-' lower(chi) 'h.stc'], output);
    
    output.data = squeeze(betas_ksdensity_amp(:,:,i));
    mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/dynamicModels/eightbins/avg_ksdensity_amp-' lower(chi) 'h.stc'], output);

    output.data = squeeze(betas_ksdensity_cf(:,:,i));
    mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/dynamicModels/eightbins/avg_ksdensity_cf-' lower(chi) 'h.stc'], output);

    output.data = squeeze(betas_normfit_sd(:,:,i));
    mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/dynamicModels/eightbins/avg_normfit_sd-' lower(chi) 'h.stc'], output);

    output.data = squeeze(betas_normfit_cf(:,:,i));
    mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/dynamicModels/eightbins/avg_normfit_cf-' lower(chi) 'h.stc'], output);
end
clear