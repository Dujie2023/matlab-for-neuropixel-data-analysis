%% load DATASTRC
np_data_pathway =  'G:\NPrecording_XY\b107a08_NP05_imec1';
filename = 'b107a08_NP05_imec1Data.mat';
load(fullfile(np_data_pathway, filename));

SpikeTimes = DataStrc.SpikeTimes;
SpikeClus =  DataStrc.SpikeClus;
UsedClusIDs = DataStrc.UsedClusIDs;
TaskTriggerStartTime = DataStrc.TaskTriggerStartTime;
behavResults = DataStrc.behavResults;

%% plot psychophysical curve
TriChoice = behavResults.Action_choice;
stim = behavResults.Stim_toneFreq;
BlockType = behavResults.BlockType;
TriType = behavResults.Trial_Type;
plot_psychophysics(np_data_pathway, TriChoice, stim, BlockType, TriType);

%%  RL model data
[P_choiceLeft,nll,delta,V_R,V_L,P_bound_low,P_bound_high, inds_correct_used,inds_use] = ...
    BSRL_djforXYdata(np_data_pathway, behavResults);

Trileft = TriChoice == 1;
TriLeft_use = Trileft(inds_use);
TriCorr_use = inds_correct_used;
blocktype_use = BlockType(inds_use);

save(fullfile(np_data_pathway,'beh_extract.mat'),'TriChoice','TriType','BlockType',...
    'stim','inds_use', 'TriLeft_use', 'TriCorr_use', 'blocktype_use');
%% extract np data
%assign spikes to clusters
spike_time_incluster = cell(length(UsedClusIDs),1);
for i = 1:length(UsedClusIDs)
    spike_time_incluster{i} = SpikeTimes(SpikeClus == UsedClusIDs(i));
end

%aline time
Tone_onset_time = behavResults.Time_stimOnset;
Answer_time = behavResults.Time_answer;
[Time_alined_outcome, Time_alined_tone,aline_time_window, aline_time_window_outcome,...
    t_after_tone_onset, t_before_tone_onset, t_after_answer, t_before_answer] = ... 
    alinetime(Tone_onset_time, Answer_time,inds_use,TaskTriggerStartTime);

%bin setting
bin_length = 10;%ms
use_ntrials = sum(inds_use);
bin_stim_onset = ceil(t_before_tone_onset/bin_length);
bin_answer = ceil(t_before_answer/bin_length);

num_bin = floor(aline_time_window/bin_length);
edges = 1:bin_length:aline_time_window;
num_bin_outcome = floor(aline_time_window_outcome/bin_length);
edges_outcome = 1:bin_length:aline_time_window_outcome;

nclusters = length(UsedClusIDs);
Spike_firingrate_bin = zeros(nclusters,use_ntrials,num_bin);% extracted np data, alined to tone onset
Spike_firingrate_bin_outcome = zeros(nclusters,use_ntrials,num_bin_outcome);% extracted np data, alined to answer
spike_time_alined = cell(use_ntrials, nclusters);
spike_time_alinedoutcome = cell(use_ntrials, nclusters);

tri_idx = find(inds_use);
Answertime_omiss = Answer_time(inds_use);
Tone_onsetT_omiss = Tone_onset_time(inds_use);
for i = 1:nclusters
     for t = 1:use_ntrials
         
        if tri_idx(t) ~= length(inds_use)
        % task time
          temp = spike_time_incluster{i}((TaskTriggerStartTime(tri_idx(t))...
            <=spike_time_incluster{i})&(spike_time_incluster{i}<=TaskTriggerStartTime(tri_idx(t)+1)));%unit:s
        else
        % for the last trial
          temp = spike_time_incluster{i}((TaskTriggerStartTime(ntrials)...
            <=spike_time_incluster{i})&(spike_time_incluster{i}<=TaskTriggerStartTime(ntrials)+10));%unit:s
        end
        
        temp_ms = 1000*(temp-TaskTriggerStartTime(tri_idx(t)));%ms from time start
        
        % aline time
        spike_time_alined{t,i} = (temp_ms((Time_alined_tone(t,1)<=temp_ms)...
            &(temp_ms<=Time_alined_tone(t,2))))-double(Tone_onsetT_omiss(t))+double(t_before_tone_onset);%ms

        %assign spikes to bins
        Spike_firingrate_bin(i,t,:) =  1000*histcounts(spike_time_alined{t,i},...
                                          edges)/bin_length;%unit: Hz
        
        % aline time to outcome
        spike_time_alinedoutcome{t,i} = (temp_ms((Time_alined_outcome(t,1)<=temp_ms)...
                &(temp_ms<=Time_alined_outcome(t,2))))-double(Answertime_omiss(t))+double(t_before_answer);%ms

        %assign spikes to bins
            Spike_firingrate_bin_outcome(i,t,:) =  1000*histcounts(spike_time_alinedoutcome{t,i},...
                                               edges_outcome)/bin_length;%unit: Hz       
     end
     
    fprintf('Cluster %d data extracted.\n', i);
    
end
save(fullfile(np_data_pathway, 'np_extracted.mat'), 'Spike_firingrate_bin','Spike_firingrate_bin_outcome', 'spike_time_alined',...
    'spike_time_alinedoutcome','bin_stim_onset','bin_answer', 't_before_tone_onset', 't_after_tone_onset', ...
    't_before_answer','t_after_answer','bin_length', 'spike_time_incluster');

%% corr analysis
np2tone = Spike_firingrate_bin;
np2answer = Spike_firingrate_bin_outcome;
BlockType = behavResults.BlockType;

modelvar = [delta; P_bound_low(2:end); P_bound_low(1:end-1); V_L; P_choiceLeft;...
            P_bound_high(2:end); P_bound_high(1:end -1)];
varname = {'Prediction error (delta)', 'P boundary low(t+1)', 'P boundary low(t)', 'Left Value', 'P choiceleft',...
            'P boundary high(t+1)', 'P bounary high(t)'};

R_allvar_bftone = zeros(size(modelvar,1),size(np2tone,1));
R_allvar_aftone = zeros(size(modelvar,1),size(np2tone,1));
R_allvar_bfanswer = zeros(size(modelvar,1),size(np2tone,1));
R_allvar_afanswer = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_bftone = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_aftone = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_bfanswer = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_afanswer = zeros(size(modelvar,1),size(np2tone,1));

save_folder = fullfile(np_data_pathway, 'corr');
if ~exist(save_folder)
    mkdir(save_folder);
end

for ivar = 1:size(modelvar,1)
    imodelvar = modelvar(ivar,:);
    ivarname = varname{ivar};
    
    R_bftone = [];
    R_aftone = [];
    R_bfanswer = [];
    R_afanswer = [];
    for iunit = 1:size(np2tone,1) 

        %before tone onset bin
        temp_bftone = mean(squeeze(np2tone(iunit,:,1:bin_stim_onset-1)),2);
        [R_bftone(iunit), P_bftone(iunit)]= corr(temp_bftone, imodelvar','Type','Pearson');

        %after tone onset
        temp_aftone = mean(squeeze(np2tone(iunit,:,bin_stim_onset:200)),2);
        [R_aftone(iunit), P_aftone(iunit)]= corr(temp_aftone, imodelvar','Type','Pearson');

        %before answer
        temp_bfanswer = mean(squeeze(np2answer(iunit,:,1:bin_answer-1)),2);
        [R_bfanswer(iunit), P_bfanswer(iunit)]= corr(temp_bfanswer, imodelvar','Type','Pearson');

        %after answer
        temp_afanswer = mean(squeeze(np2answer(iunit,:,bin_answer:end)),2);
        [R_afanswer(iunit), P_afanswer(iunit)]= corr(temp_afanswer, imodelvar','Type','Pearson');

    end
    
    R_allvar_bftone(ivar,:) = R_bftone;
    R_allvar_aftone(ivar,:) = R_aftone;
    R_allvar_bfanswer(ivar,:) = R_bfanswer;
    R_allvar_afanswer(ivar,:) = R_afanswer;
    P_allvar_bftone(ivar,:) = P_bftone;
    P_allvar_aftone(ivar,:) = P_aftone;
    P_allvar_bfanswer(ivar,:) = P_bfanswer;
    P_allvar_afanswer(ivar,:) = P_afanswer;
    
    figure;
    scatter(R_bfanswer, R_afanswer);
    xlabel('Coefficient with neural activity before answer');
    ylabel('Coefficient with neural activity after answer');
    title(['Pearson coefficient with ',ivarname]);
    saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_answer.png']);
    saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_answer.fig']);
    close;

    figure;
    scatter(R_bftone, R_aftone);
    xlabel('Coefficient with neural activity before tone');
    ylabel('Coefficient with neural activity after tone');
    title(['Pearson coefficient with ',ivarname]);
    saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_tone.png']);
    saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_tone.fig']);
    close;

    nunit = size(np2answer,1);
    per01(ivar,1:4) = [sum(P_bftone<0.01) sum(P_aftone<0.01) sum(P_bfanswer<0.01) sum(P_afanswer<0.01)]/nunit;
    per05(ivar,1:4) = [sum(P_bftone<0.05) sum(P_aftone<0.05) sum(P_bfanswer<0.05) sum(P_afanswer<0.05)]/nunit;

    % plot example neuron
    y6 = double(BlockType(inds_use));
    y6(find(y6==0))=-1;
    y6 = 1.5*y6';
    [~,example_idx] = sort(P_afanswer);
    for ie = 1:5
    figure;
    example_num = example_idx(ie);
    sgtitle(['unit ' num2str(example_num)]);
    set(gcf, 'Position', [100 100 1600 700]);
    subplot(1,2,1)
    x = 1:length(imodelvar);
    y1 = smooth(zscore(imodelvar),10);
%     y2 = smooth(zscore(mean(squeeze(np2tone(example_num,:,1:bin_stim_onset-1)),2)),10);
%     y3 = smooth(zscore(mean(squeeze(np2tone(example_num,:,bin_stim_onset:200)),2)),10);
%     y4 = smooth(zscore(mean(squeeze(np2answer(example_num,:,1:bin_answer-1)),2)),10);
    y5 = smooth(zscore(mean(squeeze(np2answer(example_num,:,bin_answer:end)),2)),10);
    plot(x,y1,'k','LineWidth',1);
    hold on;
%     plot(x,y2,'y','LineWidth',1.5);
%     plot(x,y3,'g','LineWidth',1.5);
%     plot(x,y4,'c','LineWidth',1.5);
    plot(x,y5,'b','LineWidth',1);
    plot(x,y6,'r','LineWidth',2);
    xlabel('# Trial','FontSize',12);
    ylabel('Z score','FontSize',12);
    legend('Model','Neural activity afanswer','Block','FontSize',12, 'Location', 'best');
    legend('boxoff');
    hold off;
    
    subplot(1,2,2);
    scatter(y1,y5);
    xlabel(['z-scored ',ivarname],'FontSize',12);
    ylabel('z-scored firing rate','FontSize',12);
    saveas(gcf,[np_data_pathway, '\corr\',ivarname, ' example unit ', num2str(example_num),'.png']);
    saveas(gcf,[np_data_pathway, '\corr\',ivarname, ' example unit ', num2str(example_num),'.png']);
    close;
    end
    
    %small time window
    ave_nbin =10;
    for iunit = 1:size(np2answer,1) 
        for jbin = 1:size(np2tone,3)/ave_nbin
            temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin-ave_nbin+1:jbin*ave_nbin), 3));
            [R2tone(iunit,jbin), P2tone(iunit,jbin)]= corr(temp', imodelvar','Type','Pearson');        
        end
        for jbin = 1:size(np2answer,3)/ave_nbin
            temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin-ave_nbin+1:jbin*ave_nbin), 3));
            [R2answer(iunit,jbin), P2answer(iunit,jbin)]= corr(temp', imodelvar','Type','Pearson');        
        end
    end

    % R2tone(isnan(R2tone(:))) = 0;
    % R_mean = mean(R2tone,2);
    % [~, sort_idx] = sort(R_mean);
    % for i = 1:length(R_mean)
    %     R2tone_sorted(i,:) = R2tone(sort_idx(i),:);
    % end

    R2tone_sorted = [];
    R2_answer_sorted = [];
    R2answer(isnan(R2answer(:))) = 0;
    R_mean = mean(R2answer,2);
    [~, sort_idx] = sort(R_mean);
    for i = 1:length(R_mean)
        R2answer_sorted(i,:) = R2answer(sort_idx(i),:);
        R2tone_sorted(i,:) = R2tone(sort_idx(i),:);
    end

    figure; 
    set(gcf, 'Position', [100 100 1400 700]);
    subplot(1,2,2);
    imagesc(R2answer_sorted);
    ylabel('Unit Number','FontSize',12);
    xlabel('Time (ms)','FontSize',12);
    colorbar;
    temp1 = caxis;
    set(gca, 'Xtick', [0:1000/(bin_length*ave_nbin):size(np2answer,3)/ave_nbin], 'XTicklabel', ...
        [0:1000:size(np2answer,3)*bin_length]);
    title(['Pearson coefficient with ', ivarname,' 2answer'], 'FontSize', 14);
    %answer time
    hold on;
    xanswer = bin_answer*ones(1,size(np2answer,1))/ave_nbin;
    yanswer = [1:size(np2answer,1)];
    line(xanswer, yanswer, 'LineWidth', 1.5,'LineStyle', '--', 'Color', 'k');
    hold off;

    subplot(1,2,1) 
    imagesc(R2tone_sorted(:,1:2000/(bin_length*ave_nbin)));
    caxis(temp1);
    ylabel('Unit Number','FontSize',12);
    xlabel('Time (ms)','FontSize',12);
    set(gca, 'Xtick', [0:1000/(bin_length*ave_nbin):2000/(bin_length*ave_nbin)], 'XTicklabel', ...
        [0:1000:2000]);
    title(['Pearson coefficient with ', ivarname,' 2tone'],'FontSize', 14);
    %tone time
    ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin, ...
          [size(R2tone, 1) 0 0  size(R2tone, 1)], 'w', 'edgecolor', 'none');
    alpha(ph,.4);
    filename = fullfile(np_data_pathway, 'corr', ivarname);
    saveas(gcf,[filename  '.png']);
    saveas(gcf,[filename  '.fig']);

end
save(fullfile(np_data_pathway, 'corr','corr.mat'), 'per01', 'per05','varname',...
    'R_allvar_bftone','R_allvar_aftone','R_allvar_bfanswer','R_allvar_afanswer',...
     'P_allvar_bftone','P_allvar_aftone','P_allvar_bfanswer','P_allvar_afanswer');
 
%% AUC analysis
save_folder = fullfile(np_data_pathway, 'auc');
if ~exist(save_folder)
    mkdir(save_folder);
end

blocktype_use = BlockType(inds_use);
auc_block_bftone = nan(size(np2tone,1),1);
auc_block_aftone = nan(size(np2tone,1),1);
auc_block_bfanswer = nan(size(np2answer,1),1);
auc_block_afanswer = nan(size(np2answer,1),1);

for icluster = 1:size(np2tone,1)
    neuroactivity_bftone = mean(np2tone(icluster,:,1:bin_stim_onset-1), 3);
    [auc_block_bftone(icluster), auc_block_bftone1(icluster)] = cal_auc2(neuroactivity_bftone, blocktype_use);
    
    neuroactivity_aftone = mean(np2tone(icluster,:,bin_stim_onset:200), 3);    
    [auc_block_aftone(icluster), auc_block_aftone1(icluster) ]= cal_auc2(neuroactivity_aftone, blocktype_use);
    
    neuroactivity_bfanswer = mean(np2answer(icluster,:,1:bin_answer-1), 3);
    [auc_block_bfanswer(icluster), auc_block_bfanswer1(icluster)] = cal_auc2(neuroactivity_bfanswer, blocktype_use);
    
    neuroactivity_afanswer = mean(np2answer(icluster,:,bin_answer:end), 3);    
    [auc_block_afanswer(icluster), auc_block_afanswer1(icluster)] = cal_auc2(neuroactivity_afanswer, blocktype_use);

end

figure;
scatter(auc_block_bftone,  auc_block_aftone);
xlabel('Before Tone Block AUC');
ylabel('After Tone Block AUC');
saveas(gcf, fullfile(save_folder,'bftone_vs_aftone_auc.png'));
saveas(gcf, fullfile(save_folder,'bftone_vs_aftone_auc.fig'));

figure;
scatter(auc_block_bfanswer,  auc_block_afanswer);
xlabel('Before answer Block AUC');
ylabel('After answer Block AUC');
saveas(gcf, fullfile(save_folder,'bfanswer_vs_afafanswer_auc.png'));
saveas(gcf, fullfile(save_folder,'bfanswer_vs_afafanswer_auc.fig'));

%small time window
auc2tone = [];
auc2answer = [];
ave_nbin_auc = 50;
for iunit = 1:size(np2answer,1) 
for jbin = 1:size(np2tone,3)/ave_nbin_auc
    temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin_auc-ave_nbin_auc+1:jbin*ave_nbin_auc), 3));
    [auc2tone(iunit,jbin), ~]= cal_auc2(temp, blocktype_use);        
    end
    for jbin = 1:size(np2answer,3)/ave_nbin_auc
        temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin_auc-ave_nbin_auc+1:jbin*ave_nbin_auc), 3));
        [auc2answer(iunit,jbin),~]= cal_auc2(temp, blocktype_use);        
    end
end

auc2tone_sorted = [];
auc2answer_sorted = [];
auc2answer(isnan(auc2answer(:))) = 0;
auc_mean = mean(auc2answer,2);
[~, sort_idx] = sort(auc_mean);
for i = 1:length(auc_mean)
    auc2answer_sorted(i,:) = auc2answer(sort_idx(i),:);
    auc2tone_sorted(i,:) = auc2tone(sort_idx(i),:);
end

figure; 
set(gcf, 'Position', [100 100 1400 700]);
subplot(1,2,2);
imagesc(auc2answer_sorted);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
colorbar;
temp1 = caxis;
set(gca, 'Xtick', [0:1000/(bin_length*ave_nbin_auc):size(np2answer,3)/ave_nbin_auc], 'XTicklabel', ...
   [0:1000:size(np2answer,3)*bin_length]);
title(['Block AUC ', ' 2answer'], 'FontSize', 14);
%answer time
hold on;
xanswer = bin_answer*ones(1,size(np2answer,1))/ave_nbin_auc;
yanswer = [1:size(np2answer,1)];
line(xanswer, yanswer, 'LineWidth', 1.5,'LineStyle', '--', 'Color', 'k');
hold off;

subplot(1,2,1) 
imagesc(auc2tone_sorted(:,1:2000/(bin_length*ave_nbin_auc)));
caxis(temp1);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
set(gca, 'Xtick', [0:1000/(bin_length*ave_nbin_auc):2000/(bin_length*ave_nbin_auc)], 'XTicklabel', ...
    [0:1000:2000]);
title(['Block AUC ',' 2tone'],'FontSize', 14);
%tone time
ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin_auc, ...
      [size(auc2tone, 1) 0 0  size(auc2tone, 1)], 'w', 'edgecolor', 'none');
alpha(ph,.4);
filename = fullfile(np_data_pathway, 'auc','allunit');
saveas(gcf,[filename  '.png']);
saveas(gcf,[filename  '.fig']);

% figure;
% scatter(2*(auc_block_bftone-0.5),  2*(auc_block_aftone-0.5));
% xlabel('Baseline Block Selectivity Index');
% ylabel('After Baseline Block Selectivity Index');
% saveas(gcf, 'baseline_vs_afbaseline_SI.png');
save(fullfile(np_data_pathway, 'auc','auc.mat'), 'auc_block_bftone', 'auc_block_aftone',...
'auc_block_bfanswer', 'auc_block_afanswer','auc_block_bftone1', 'auc_block_aftone1',...
'auc_block_bfanswer1', 'auc_block_afanswer1');

