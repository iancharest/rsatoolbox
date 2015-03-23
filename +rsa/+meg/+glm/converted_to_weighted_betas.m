path = '/imaging/iz01/NeuroLex/glm/sixteenbins/';
load('weights_16bins_256numchans_1-10kHz.mat');
chi = 'R';
nbins = 16;

scale = [-100; -50; -25; -12.5; -6.2; -3.1; -1.5; -0.7; 0.7; 1.5; 3.1; 6.2; 12.5; 25; 50; 100]; 
% models = {'bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'bin9', 'bin10', 'bin11', 'bin12', 'bin13', 'bin14', 'bin15', 'bin16'};
k=1;bin_weight=[];
for i=1:16:256
    bin_weight(k) = mean(weights(i:i+15));
    k = k+1;
    disp(num2str(i));
end
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

for i = 1:length(subjectNames)
    disp(['subject: ' num2str(i)]);
    betas = mne_read_stc_file1([path 'avg_glm_beta-' lower(chi) 'h.stc']);
    bins = mne_read_stc_file1([path 'avg_glm_scale-' lower(chi) 'h.stc']);

    scaled = betas;
    for j = 1:nbins
        scaled.data((round(bins.data*10)/10)==scale(j)) = betas.data((round(bins.data*10)/10)==scale(j)) .* bin_weight(j);
        if isempty(betas.data((round(bins.data*10)/10)==scale(j)))
            disp([' bin ' num2str(j) ' not found!']);
        end
        if find(betas.data(bins.data==scale(j))<0)
            disp([' Negative values in bin ' num2str(j)]);
        end
    end    
     mne_write_stc_file1([path 'avg_weighted_glm_beta-' lower(chi) 'h.stc'], scaled);
end