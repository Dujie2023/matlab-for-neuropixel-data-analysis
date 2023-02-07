function [cluster_feature, spike_time_incluster] = cal_filter_features(np_data_pathway, ksfolder)
%np_data_pathway =  'G:\NPrecording_XY\Cat_A20220101_b107a07_NPSess06_g0_imec1';

%load data
spike_times = readNPY(fullfile(ksfolder, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(ksfolder, 'spike_clusters.npy'));
cluster_info_path = fullfile(ksfolder,'cluster_info.csv');
%assign spikes to clusters
clusters = unique(spike_clusters);
spike_time_incluster = cell(length(clusters),1);
for i = 1:length(clusters)
    spike_time_incluster{i} = spike_times(spike_clusters == clusters(i));
end

cluster_info = readtable(cluster_info_path );
cluster_id = cluster_info.id;

load(fullfile(ksfolder,'WaveformDatas.mat'),'UnitFeatures');
%% calculating metrics
%kslabel
cluster_feature.kslabel = cellfun(@(x) strcmp(x,'good'), cluster_info.KSLabel);%logic mat(nclusters,1)

%firing rate
%cluster_feature.firing_rate = cluster_info.fr;%matrix (nclusters,1)
cluster_feature.fr = cal_fr(np_data_pathway, spike_time_incluster);% fr in task time, unit: Hz

% presence ratio
nbins = 20;
cluster_feature.presence_ratio = cal_presence_ratio(np_data_pathway, spike_time_incluster,nbins);%matrix (nclusters,1)

% ISI violations
minISI = 1*30;
refDur = 1.5*30;
[cluster_feature.ISI_violations, numviolations] =  cal_ISIViolations(spike_time_incluster, minISI, refDur);

%unitfeature(num_cluster,{wave_amplitude, trough2peakT, repolarT,
%                       prepeak2pospeakratio, wave_spread, snr, troughpeakidx});

%waveform_ampltitude
cluster_feature.wave_ampltitude = cell2mat(UnitFeatures(:,1));

%trough2peakT
cluster_feature.trough2peakT = cell2mat(UnitFeatures(:,2))/30;%unit: ms

%repolarization time
cluster_feature.repolarT = cell2mat(UnitFeatures(:,3))/30;%unit: ms

%prepeak2pospeakratio
cluster_feature.prepeak2pospeakratio = cell2mat(UnitFeatures(:,4));

%waveform spread
cluster_feature.wave_spread = cell2mat(UnitFeatures(:,5));%unit: channel

%signal-to-noise ratio
cluster_feature.snr = cell2mat(UnitFeatures(:,6));

save(fullfile(ksfolder, 'filter_features,mat'), 'cluster_feature', 'spike_time_incluster');
end

function [fr] = cal_fr(np_path, spike_time_incluster)
%the presence ratio is defined as the fraction of blocks that include one or more spikes from a particular unit. 

load(fullfile(np_path, 'TriggerDatas.mat'));

%task time
time_after_finaltrigger = 20;
samplingrate = 30000;
task_start_time = TriggerEvents(1,1);
task_end_time = TriggerEvents(end,1) + samplingrate*time_after_finaltrigger;
T = task_end_time - task_start_time;

fr = cellfun(@(x) sum((task_start_time<=x)&(x<=task_end_time))*samplingrate /T, spike_time_incluster);

end

function [presence_ratio] = cal_presence_ratio(np_path, spike_time_incluster,nbins)
%the presence ratio is defined as the fraction of blocks that include one or more spikes from a particular unit. 

load(fullfile(np_path, 'TriggerDatas.mat'));

%task time
time_after_finaltrigger = 20;
samplingrate = 30000;
task_start_time = TriggerEvents(1,1);
task_end_time = TriggerEvents(end,1) + samplingrate*time_after_finaltrigger;

edges = round(linspace(task_start_time, task_end_time, nbins+1));

%calculate presence presence ratio
presence_ratio = nan(length(spike_time_incluster),1);
for ncluster = 1:length(spike_time_incluster)
    temp = histcounts(spike_time_incluster{ncluster},edges);% matrix(nbin,1), count spike numers in each bin
    presence_ratio(ncluster) = (sum(temp>=1))/nbins;
end

end

function [fpRate, numViolations] = cal_ISIViolations(spike_time_incluster, minISI, refDur)

fpRate = nan(length(spike_time_incluster),1);
numViolations = nan(length(spike_time_incluster),1);

for ncluster = 1:length(spike_time_incluster)
        [fpRate(ncluster),numViolations(ncluster)] = ISIViolations...
            (double(spike_time_incluster{ncluster}), minISI, refDur);% matrix(nbin,1), count spike numers in each bin;
end

end

function [fpRate, numViolations] = ISIViolations(spikeTrain, minISI, refDur)
% computes an estimated false positive rate of the spikes in the given
% spike train. You probably don't want this to be greater than a few
% percent of the spikes. 
%
% - minISI [sec] = minimum *possible* ISI (based on your cutoff window); likely 1ms
% or so. It's important that you have this number right.
% - refDur [sec] = estimate of the duration of the refractory period; suggested value = 0.002.
% It's also important that you have this number right, but there's no way
% to know... also, if it's a fast spiking cell, you could use a smaller
% value than for regular spiking.

% The false positive value should less than 0.1, Ref from 
% https://doi.org/10.1038/s41586-020-03166-8

totalRate = length(spikeTrain)/spikeTrain(end);
numViolations = sum(diff(spikeTrain) <= refDur);

% Spikes occurring under the minimum ISI are by definition false positives;
% their total rate relative to the rate of the neuron overall gives the
% estimated false positive rate.
% NOTE: This does not use the equation from Dan Hill's paper (see below)
% but instead uses the equation from his UltraMegaSort
violationTime = 2*length(spikeTrain)*(refDur-minISI); % total time available for violations - 
                                                    % there is an available
                                                    % window of length
                                                    % (refDur-minISI) after
                                                    % each spike.
violationRate = numViolations/violationTime;
fpRate = violationRate/totalRate;

if fpRate>1
    % it is nonsense to have a rate >1, however having such large rates
    % does tell you something interesting, namely that the assumptions of
    % this analysis are failing!
    fpRate = 1; 
end


% Here is with Dan Hill's equation (J. Neurosci, 2012)
%  NOTE: actually, this method will return imaginary results for large-ish
%  numbers of violations. 
% tauR = refDur;
% tauC = minISI;
% T = spikeTrain(end);
% N = length(spikeTrain);
% R = sum(diff(spikeTrain) <= tauR);
% 
% k = 2*(tauR-tauC)*N^2;
% rts = roots([k -k R*T]);
% 
% fpRate = min(rts);
% numViolations = R;
end