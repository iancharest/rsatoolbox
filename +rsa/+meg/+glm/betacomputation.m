toolboxRoot = '/imaging/iz01/test/toolbox/devel/toolbox'; addpath(genpath(toolboxRoot));
% loading response
path = '/imaging/ls02/NeuroLex/cochleagram_HD/RDMs/';
chi = 'R';
load([path 'averaged_searchlightRDMs_masked_' lower(chi) 'h.mat']);

% loading predictors
load('Models-GTF1to30ms.mat');
% load('VoiceModel.mat');
% load('placeModel.mat');
% load('MannerModel.mat');

% load('Low-Half.mat');
% load('High-Half.mat');
load('Low1.mat');
load('Low2.mat');

load('Mid1.mat');
load('Mid2.mat');
load('Mid3.mat');
load('Mid4.mat');

load('High1.mat');
load('High2.mat');

betas = zeros(10242,22); 
% scale for coloring
% scale = [15; 5; -5; -15]; % cochlear, voice, place, manner
% models = {'cochlear', 'voice', 'place', 'manner'};
% scale = [-15; 15]; % low, high
% models = {'low', 'high'};
% scale = {-1; 1};
scale = [-15; -7.5; -3.7; -1.8; 1.8; 3.7; 7.5; 15]; 
models = {'lo1', 'lo2', 'mid1', 'mid2', 'mid3', 'mid4', 'high1', 'high2'};


% applying glm
vertices = fieldnames(averageSubjectRDMs.(chi));

for v =1:length(vertices)
    thisVertex = vertices{v}; 
    [temp vertex_num] = strtok(thisVertex,'_'); vertex_num = str2num(strtok(vertex_num,'_'));
    nTimePoints = fieldnames(averageSubjectRDMs.(chi).(thisVertex));
    
    for t = 1:length(nTimePoints) 
        thisTime = nTimePoints{t};

%         [b{v,t},dev{v,t},stats{v,t}] = glmfit([vectorizeRDM(low_half.RDM); vectorizeRDM(high_half.RDM)]',...
%             averageSubjectRDMs.(chi).(thisVertex).(thisTime).RDM','normal');
%         
        [b{v,t},dev{v,t},stats{v,t}] = glmfit([vectorizeRDM(low1.RDM);...
            vectorizeRDM(low2.RDM); vectorizeRDM(mid1.RDM); ...
            vectorizeRDM(mid2.RDM); vectorizeRDM(mid3.RDM); ...
            vectorizeRDM(mid4.RDM); vectorizeRDM(high1.RDM); ...
            vectorizeRDM(high2.RDM)]',...
            averageSubjectRDMs.(chi).(thisVertex).(thisTime).RDM','normal');
        
        which_model = find(b{v,t}==max(b{v,t}(2:length(b{v,t}),:))) - 1; % 1's are added by glm as predictor
        
%         betas(vertex_num,t) = scale(which_model) * max(b{v,t}(2:length(b{v,t}),:));
        betas(vertex_num,t) = scale(which_model);
       
        disp([thisVertex ' ' thisTime ' ' models{which_model}]);
        
    end    
end

% output settings

path = '/imaging/iz01/NeuroLex/cochleagram_HD/Results/';

output = mne_read_stc_file1([path 'lexPro_tonotopy_1to30ms_HD_Gamma30ms_significant_vertex-lh.stc']);
output.tmin = -0.05;
output.tstep = 0.01;
output.data = betas;

mne_write_stc_file1(['/imaging/iz01/NeuroLex/glm/LoHiFreq/eightbins-' lower(chi) 'h.stc'], output);

save('glmResults','output','b','scale','models','dev','stats');