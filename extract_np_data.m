function [Spike_firingrate_bin,Spike_firingrate_bin_outcome, spike_time_alined, spike_time_alinedoutcome, bin_stim_onset, bin_answer, bin_length] = ...
            extract_np_data(np_data_pathway, behfile, alined2outcome)
%% load NP data
%np_data_pathway =  'G:\NPrecording_XY\Cat_A20220101_b107a07_NPSess06_g0_imec1';
ksfolder = fullfile(np_data_pathway,'ks2_5');

spike_times = readNPY(fullfile(ksfolder, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(ksfolder, 'spike_clusters.npy'));
cgsFile = fullfile(ksfolder,'cluster_info.csv');
cluster_info = readtable(cgsFile);
UsedClus_IDs = cluster_info.id;

load(fullfile(np_data_pathway, 'TriggerDatas.mat'));

%assign spikes to clusters
clusters = unique(spike_clusters);
spike_time_incluster = cell(length(clusters),1);
for i = 1:length(clusters)
    spike_time_incluster{i} = spike_times(spike_clusters == clusters(i));
end

%% load behavior data
%behfile = 'A20220101_b107a07_NPSess06_2afc.mat';
load(fullfile(np_data_pathway, behfile));
ntrials = size(SessionResults,2); % number of trials

%if num_trials in behavior is different from num_trials in np data
if size(TriggerEvents,1) ~= ntrials
    fprintf('Trial number in beh data and Trigger data is different.\n');
    pause;
end

Tone_onset_time = cellfun(@(x) x.Time_stimOnset, SessionResults);
TriChoice = cellfun(@(x) x.Action_choice, SessionResults);
TriType = cellfun(@(x) x.Trial_Type, SessionResults);
Miss_Ind = cellfun(@(x) x.Action_choice == 2, SessionResults);
Answer_time = cellfun(@(x) x.Time_answer, SessionResults);
Reward_time = cellfun(@(x) x.Time_reward, SessionResults);
BlockType = cellfun(@(x) x.BlockType, SessionResults);
Stim_toneFreq = cellfun(@(x) x.Stim_toneFreq, SessionResults);

Animal_Name_Setting = SessionSettings{1}.animalName;
Answer_Period = SessionSettings{1}.answerPeriod;
stimDuration = SessionSettings{1,1}.stimDuration;
FixedLeft_tones = SessionSettings{1, 1}.FixedLeft_tones;
FixedRight_tones = SessionSettings{1, 1}.FixedRight_tones;
Middleshift_tones = SessionSettings{1, 1}.Middleshift_tones;

inds_use = ~Miss_Ind;
Trileft = TriChoice == 1;
TriLeft_use = Trileft(inds_use);
blocktype_use = BlockType(inds_use);

inds_correct = Reward_time;
inds_correct_used = inds_correct(inds_use);
TriCorr_use = inds_correct_used;

MinOnsetTime = min(Tone_onset_time(inds_use));

save(fullfile(np_data_pathway,'beh_extract.mat'),'ntrials','Tone_onset_time','TriChoice','TriType','Miss_Ind',...
    'Reward_time', 'Answer_time','BlockType','Stim_toneFreq','Animal_Name_Setting',...
    'Answer_Period','stimDuration','FixedLeft_tones','FixedRight_tones',...
    'Middleshift_tones','MinOnsetTime','inds_correct_used',...
    'inds_use', 'TriLeft_use', 'TriCorr_use', 'blocktype_use');

%% aline time to tone
t_before_tone_onset = 1000;%ms
t_after_tone_onset = 5000;%ms
if MinOnsetTime <= t_before_tone_onset
    t_before_tone_onset = double(MinOnsetTime);
end

Tone_onsetT_omiss = Tone_onset_time(inds_use);
trial_all_time = (TriggerEvents(2:ntrials+1,1)-TriggerEvents(1:ntrials,1))/30;
min_tone2endT = round(min(trial_all_time(inds_use)' - double(Tone_onsetT_omiss)));
if t_after_tone_onset >= min_tone2endT
    t_after_tone_onset = min_tone2endT;
end

for triN = 1:length(Tone_onsetT_omiss)
    Time_alined_tone(triN,:) = double(Tone_onsetT_omiss(triN)) + [-t_before_tone_onset+1 t_after_tone_onset];% time 1s before tone onset, time 6-8s after tone offset;
end
aline_time_window =  t_after_tone_onset + t_before_tone_onset;

%aline time to outcome
t_before_answer = 1000;
t_after_answer  = 4000;
Answertime_omiss = Answer_time(inds_use);
MinAnswerT = min(Answertime_omiss);
if MinAnswerT <= t_before_answer
    t_before_answer = MinAnswerT;
end
MinAnswer2endT = round(min(trial_all_time(inds_use)' - double(Answertime_omiss)));
if MinAnswer2endT <= t_after_answer
    t_after_answer = MinAnswer2endT;
end

for triN = 1:length(Answertime_omiss)
    Time_alined_outcome(triN,:) = double(Answertime_omiss(triN)) + [-t_before_answer+1 t_after_answer];
end
aline_time_window_outcome =  t_after_answer + t_before_answer;
%% assign spikes to bins
 
%bin setting
bin_length = 10;%ms
use_ntrials = sum(inds_use);
num_bin = floor(aline_time_window/bin_length);
edges = 1:bin_length:aline_time_window;

num_bin_outcome = floor(aline_time_window_outcome/bin_length);
edges_outcome = 1:bin_length:aline_time_window_outcome;

Spike_firingrate_bin = zeros(length(clusters),use_ntrials,num_bin);% extracted np data, alined to tone onset
Spike_firingrate_bin_outcome = zeros(length(clusters),use_ntrials,num_bin_outcome);% extracted np data, alined to answer
spike_time_alined = cell(use_ntrials, length(clusters));
spike_time_alinedoutcome = cell(use_ntrials, length(clusters));

bin_stim_onset = ceil(t_before_tone_onset/bin_length);
bin_answer = ceil(t_before_answer/bin_length);
tri_idx = find(inds_use);
for i = 1:length(clusters)
     for t = 1:use_ntrials
         
        if tri_idx(t) ~= ntrials
        % task time
          temp = spike_time_incluster{i}((TriggerEvents(tri_idx(t),1)...
            <=spike_time_incluster{i})&(spike_time_incluster{i}<=TriggerEvents(tri_idx(t)+1,1)));
        else
        % for the last trial
          temp = spike_time_incluster{i}((TriggerEvents(ntrials,1)...
            <=spike_time_incluster{i})&(spike_time_incluster{i}<=TriggerEvents(ntrials,1)+30000*20));
        end
        
        temp_ms = double(temp-TriggerEvents(tri_idx(t),1))/30;%ms from time start
        
        % aline time
        spike_time_alined{t,i} = (temp_ms((Time_alined_tone(t,1)<=temp_ms)...
            &(temp_ms<=Time_alined_tone(t,2))))-double(Tone_onsetT_omiss(t))+double(t_before_tone_onset);%ms

        %assign spikes to bins
        Spike_firingrate_bin(i,t,:) =  1000*histcounts(spike_time_alined{t,i},...
                                          edges)/bin_length;%unit: Hz
        
        if alined2outcome
            % aline time to outcome
            spike_time_alinedoutcome{t,i} = (temp_ms((Time_alined_outcome(t,1)<=temp_ms)...
                &(temp_ms<=Time_alined_outcome(t,2))))-double(Answertime_omiss(t))+double(t_before_answer);%ms

            %assign spikes to bins
            Spike_firingrate_bin_outcome(i,t,:) =  1000*histcounts(spike_time_alinedoutcome{t,i},...
                                                  edges_outcome)/bin_length;%unit: Hz
        end
        
     end
     
    fprintf('Cluster %d data extracted.\n', i);
    
end
save(fullfile(np_data_pathway, 'np_extracted.mat'), 'Spike_firingrate_bin','Spike_firingrate_bin_outcome', 'spike_time_alined',...
    'spike_time_alinedoutcome','bin_stim_onset','bin_answer','stimDuration', 't_before_tone_onset', 't_after_tone_onset', ...
    'MinAnswerT','t_before_answer','t_after_answer','bin_length', 'UsedClus_IDs','spike_time_incluster');


end
