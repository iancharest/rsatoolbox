toolboxRoot = '/imaging/iz01/test/toolbox/devel/toolbox'; addpath(genpath(toolboxRoot));
% loading response
path = '/imaging/ls02/NeuroLex/cochleagram_HD/RDMs/';
chi = 'L';
nbins = 16;
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

temppath = '/imaging/iz01/NeuroLex/cochleagram_HD/Results/';

output = mne_read_stc_file1([temppath 'lexPro_tonotopy_1to30ms_HD_Gamma30ms_significant_vertex-lh.stc']);
output.tmin = -0.05;
output.tstep = 0.01;

betas_glm_amp = zeros(10242,22);
betas_glm_scale = zeros(10242,22);

betas_ksdensity_amp = zeros(10242,22);
betas_ksdensity_cf = zeros(10242,22);

betas_normfit_sd = zeros(10242,22);
betas_normfit_cf = zeros(10242,22);

in_path= '/imaging/iz01/NeuroLex/glm/sixteenbins/';

for i=1:length(subjectNames)
    disp(subjectNames{i});
    
    data = mne_read_stc_file1([in_path 'glm_beta-' subjectNames{i} '-' lower(chi) 'h.stc']);
    betas_glm_amp = betas_glm_amp + data.data;
    clear data.data;
    
    data = mne_read_stc_file1([in_path 'glm_scale-' subjectNames{i} '-' lower(chi) 'h.stc']);
    betas_glm_scale = betas_glm_scale + data.data;
    clear data.data;
    
    
    data = mne_read_stc_file1([in_path 'ksdensity_amp-' subjectNames{i} '-' lower(chi) 'h.stc']);
    betas_ksdensity_amp = betas_ksdensity_amp + data.data;
    clear data.data;
    
    data = mne_read_stc_file1([in_path 'ksdensity_cf-' subjectNames{i} '-' lower(chi) 'h.stc']);
    betas_ksdensity_cf = betas_ksdensity_cf + data.data;
    logged = log(data.data);
    data.data = logged;
    data.data(data.data==-Inf) = 0;
    mne_write_stc_file1([in_path 'logged/ksdensity_cf_logged-' subjectNames{i} '-' lower(chi) 'h.stc'],data);
    clear data.data;
    
    data = mne_read_stc_file1([in_path 'normfit_sd-' subjectNames{i} '-' lower(chi) 'h.stc']);
    betas_normfit_sd = betas_normfit_sd + data.data;
    clear data.data;
    
    data = mne_read_stc_file1([in_path 'normfit_cf-' subjectNames{i} '-' lower(chi) 'h.stc']);
    betas_normfit_cf = betas_normfit_cf + data.data;
    logged = log(data.data);
    data.data = logged;
    data.data(data.data==-Inf) = 0;
    mne_write_stc_file1([in_path 'logged/normfit_cf_logged-' subjectNames{i} '-' lower(chi) 'h.stc'],data);
    clear data.data;
    
end

betas_glm_amp = betas_glm_amp/length(subjectNames);
data.data = betas_glm_amp;
mne_write_stc_file1([in_path 'avgd/glm_beta_avgd-' lower(chi) 'h.stc'],data);
clear data.data;

betas_glm_scale = betas_glm_scale/length(subjectNames);
data.data = betas_glm_scale;
mne_write_stc_file1([in_path 'avgd/glm_scale_avgd-' lower(chi) 'h.stc'],data);
clear data.data;


betas_ksdensity_amp = betas_ksdensity_amp/length(subjectNames);
data.data = betas_ksdensity_amp;
mne_write_stc_file1([in_path 'avgd/ksdensity_amp_avgd-' lower(chi) 'h.stc'],data);
clear data.data;

betas_ksdensity_cf = betas_ksdensity_cf/length(subjectNames);
data.data = betas_ksdensity_cf;
mne_write_stc_file1([in_path 'avgd/ksdensity_cf_avgd-' lower(chi) 'h.stc'],data);
data.data = log(data.data);
data.data(data.data==-Inf) = 0;
mne_write_stc_file1([in_path 'logged/ksdensity_cf_avgd_logged-' lower(chi) 'h.stc'],data);
clear data.data;


betas_normfit_sd = betas_normfit_sd/length(subjectNames);
data.data = betas_normfit_sd;
mne_write_stc_file1([in_path 'avgd/normfit_sd_avgd-' lower(chi) 'h.stc'],data);
clear data.data;

betas_normfit_cf = betas_normfit_cf/length(subjectNames);
data.data = betas_normfit_cf;
mne_write_stc_file1([in_path 'avgd/normfit_cf_avgd-' lower(chi) 'h.stc'],data);
data.data = log(data.data);
data.data(data.data==-Inf) = 0;
mne_write_stc_file1([in_path 'logged/normfit_cf_avgd_logged-' lower(chi) 'h.stc'],data);
clear data.data;
