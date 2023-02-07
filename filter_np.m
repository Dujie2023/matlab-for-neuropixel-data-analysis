function [filter_idx] = filter_np(np_data_pathway)
%input:
%output: 

np_data_pathway =  'G:\NPrecording_XY\Cat_A20220101_b107a07_NPSess06_g0_imec1';
ksfolder = fullfile(np_data_pathway,'ks2_5');

[cluster_feature, spike_time_incluster] = cal_filter_features(np_data_pathway, ksfolder);

%% defining criteria for different brain area
fr_threshold = 0.1;
ISI_threshold = 0.2;
am_threshold = 10;
prepeak2pospeakratio_threshold = 0.8;
wave_spread_threshold = 20;
snr_threshold = 1;
trough2peakT_lb = 4/30;trough2peakT_ub = 40/30;
repolarT_lb = 2/30;repolarT_ub = 30/30;

%% filtering
numcluster = size(cluster_feature.kslabel,1);
%kslabel
fl.kslabel_idx = (cluster_feature.kslabel == 1);
fl.kslabel_ratio = 1 - sum(fl.kslabel_idx)/numcluster;

%firing rate
fl.fr_idx = (cluster_feature.fr >= fr_threshold);
fl.fr_ratio =  1 - sum(fl.fr_idx)/numcluster;

% presence ratio
fl.presence_ratio_idx = (cluster_feature.presence_ratio == 1);
fl.presence_ratio_ratio = 1 - sum(fl.presence_ratio_idx)/numcluster;

% ISI violations
fl.ISI_idx =  (cluster_feature.ISI_violations <= ISI_threshold);
fl.ISI_ratio = 1 - sum(fl.ISI_idx)/numcluster;

% waveform_ampltitude
fl.am_idx = (cluster_feature.wave_ampltitude >= am_threshold);
fl.am_ratio = 1 - sum(fl.am_idx)/numcluster;

% trough2peakT
fl.trough2peakT_idx = (cluster_feature.trough2peakT >= trough2peakT_lb)&...
            (cluster_feature.trough2peakT <= trough2peakT_ub);
fl.trough2peakT_ratio = 1 - sum(fl.trough2peakT_idx)/numcluster;

% repolarization time
fl.repolarT_idx = (cluster_feature.repolarT >= repolarT_lb)&...
            (cluster_feature.repolarT <= repolarT_ub);
fl.repolarT_ratio = 1 - sum(fl.repolarT_idx)/numcluster;

% prepeak2pospeakratio
fl.prepeak2pospeakratio_idx = (cluster_feature.prepeak2pospeakratio <= prepeak2pospeakratio_threshold);
fl.prepeak2pospeakratio_ratio = 1 - sum(fl.prepeak2pospeakratio_idx)/numcluster;

% waveform spread
fl.wave_spread_idx = (cluster_feature.wave_spread <= wave_spread_threshold);
fl.wave_spread_ratio = 1 - sum(fl.wave_spread_idx)/numcluster;

% signal-to-noise ratio
fl.snr_idx = (cluster_feature.snr >= snr_threshold);
fl.snr_ratio = 1 - sum(fl.snr_idx)/numcluster;

fl.flidx = fl.kslabel_idx & fl.fr_idx & fl.presence_ratio_idx & ...
           fl.am_idx & fl.ISI_idx & fl.trough2peakT_idx & fl.repolarT_idx & ...
           fl.prepeak2pospeakratio_idx & fl.wave_spread_idx & fl.snr_idx;
fl.fl_ratio = 1 - sum(fl.flidx) / numcluster;

filter_idx = fl.flidx;

save(fullfile(np_data_pathway, 'filter.mat'), 'fl');
%% showing the process of filtering


end

