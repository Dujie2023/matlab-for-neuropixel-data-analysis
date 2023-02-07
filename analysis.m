%% extract_waveformdata
rez_path = fullfile(ksfolder,'rez2.mat');
% load(rez_path);
bin_name = A20220101_b107a07_NPSess06_g0_tcat.imec1.ap.bin;
% extract_waveformdata_dj(np_data_pathway, bin_name, rez);

%% extract np data
%np_data_pathway =  'G:\NPrecording_XY\Cat_A20220101_b107a07_NPSess06_g0_imec1';
%behfile = 'A20220101_b107a07_NPSess06_2afc.mat';
[Spike_firingrate_bin, Spike_firingrate_bin_outcome, spike_time_alined,spike_time_alinedoutcome, bin_stim_onset, bin_answer, bin_length] = ...
    extract_np_data(np_data_pathway,behfile, 1);%(ncluster,ntrial, nbins)

%% filter
load(fullfile(np_data_pathway, 'np_extracted.mat'), 'Spike_firingrate_bin');
load(fullfile(np_data_pathway, 'np_extracted.mat'), 'Spike_firingrate_bin_outcome');
load(fullfile(np_data_pathway, 'np_extracted.mat'), 'UsedClus_IDs');

% filter_idx = filter_np(np_data_pathway);
load(fullfile(np_data_pathway, 'filter.mat'));
flidx = fl.flidx;

np2tone = Spike_firingrate_bin(flidx,:,:);
np2answer = Spike_firingrate_bin_outcome(flidx,:,:);

save(fullfile(np_data_pathway, 'npdata_alined_filed.mat'), 'np2tone', 'np2answer');

%% filter (Xin Yu)
load(fullfile(np_data_pathway, 'np_extracted.mat'), 'Spike_firingrate_bin');
load(fullfile(np_data_pathway, 'np_extracted.mat'), 'Spike_firingrate_bin_outcome');
load(fullfile(np_data_pathway, 'np_extracted.mat'), 'UsedClus_IDs');

[cluster_Ab, usedAb, flidxAb] = extract_brainarea(SessAreaIndexStrc, UsedClus_IDs);
ClusIDs_filtered = UsedClus_IDs(flidxAb);
ClusAbs_filtered = cluster_Ab(flidxAb);

Abs = unique(ClusAbs_filtered);
Ab_idx = cell(1,length(Abs));
for iAb = 1:length(Abs)
    Ab_idx{iAb} = find(ClusAbs_filtered == iAb);
end

np2tone = Spike_firingrate_bin(flidxAb,:,:);
np2answer = Spike_firingrate_bin_outcome(flidxAb,:,:);
save(fullfile(np_data_pathway, 'SessNeuralData.mat'), 'np2tone', 'np2answer', 'Ab_idx',...
    'ClusIDs_filtered', 'usedAb', 'ClusAbs_filtered');

%%  RL model data
[P_choiceLeft,nll,delta,V_R,V_L,P_bound_low,P_bound_high, inds_correct_used,inds_use] = BS_RL_script_dj(np_data_pathway, behfile);

%% corr analysis
modelvar = [delta; P_bound_low(2:end); P_bound_low(1:end-1); V_L; P_choiceLeft;...
            P_bound_high(2:end); P_bound_high(1:end -1)];
varname = {'Prediction error (delta)', 'P boundary low(t+1)', 'P boundary low(t)', 'Left Value', 'P choiceleft',...
            'P boundary high(t+1)', 'P bounary high(t)'};
corr_analysis;

%% auc analysis
