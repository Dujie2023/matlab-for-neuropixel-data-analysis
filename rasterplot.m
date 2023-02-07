%% raster plot and psth
np_data_pathway = 'D:\DJ_data\catprocessdata\202206batch_mouse1_cat\Cat_20220826_mouse1_fourSites_g0_imec0\ks2_5';
load(fullfile(np_data_pathway, 'np_extracted.mat'));
np2tone = Spike_firingrate_bin;
ntrials= size(np2tone, 1);
nunit = size(np2tone, 2);
numbin = size(np2tone, 3);
color_type = {[0.118, 0.565, 1], [0.529, 0.808, 0.980], [1, 0, 0], [1, 0.39, 0.28]};

save_rasterplot_folder = fullfile(np_data_pathway, 'raster_plot');
if ~exist(save_rasterplot_folder)
    mkdir(save_rasterplot_folder);
end

%idx for different types of trials
load(fullfile(np_data_pathway,'BehavData.mat'));
inds_use = behavResults.Action_choice~=2;
blocktype_use = behavResults.BlockType(inds_use);
choice_use = behavResults.Action_choice(inds_use);
outcome = behavResults.Time_reward~=0;
inds_correct_used = outcome(inds_use);
fre_use = behavResults.Stim_toneFreq(inds_use);
fre_type = sort(unique(fre_use));
fre_corr_block_idx = cell(length(fre_type), 4);

for i = 1 : length(fre_type)
    fre_corr_block_idx{i,1} = fre_use == fre_type(i) & blocktype_use == 0 & inds_correct_used == 1;
    fre_corr_block_idx{i,2} = fre_use == fre_type(i) & blocktype_use == 0 & inds_correct_used == 0;
    fre_corr_block_idx{i,3} = fre_use == fre_type(i) & blocktype_use == 1 & inds_correct_used == 1;
    fre_corr_block_idx{i,4} = fre_use == fre_type(i) & blocktype_use == 1 & inds_correct_used == 0;
end

for iunit = 1%:nunit
    figure;
%    savename = ['Cluster ',num2str(Used_id(iunit)), ' plot'];
%    sgtitle(savename);
    set(gcf, 'Position', [100 100 1800 1100]);
    for ifre = 1 : length(fre_type)
        for j = 1:4
            trials_unit_spiketime = np_data_cell(fre_corr_block_idx{ifre,j},iunit);
            if length(trials_unit_spiketime) >= 3
                subplot(5,length(fre_type), ifre + (j-1) * length(fre_type))
                for t = 1:length(trials_unit_spiketime)
                    t_unit_spike = trials_unit_spiketime{t}; % Spike timings in the t_th trial for unit_idx
                    nspikes = numel(t_unit_spike); % number of spikes
                    for ii = 1:nspikes % for every spike
                        pl = line([t_unit_spike(ii) t_unit_spike(ii)],[t-1 t] );
                        pl.Color = 'Black';% draw a black vertical line of length 1 at time t(x) and at trial t (y)
                    end
                end
                ylim([0,length(trials_unit_spiketime)]);
                xlim([0,bin_length*size(np2tone,3)]);
                if j == 1
                    title(num2str(fre_type(ifre)));
                end
                if ifre == 1
                    if j == 1
                        ylabel('Blocktype 0 | Correct', 'Color', color_type{j});
                    elseif j == 2
                        ylabel('Blocktype 0 | Error', 'Color', color_type{j});
                    elseif j == 3
                        ylabel('Blocktype 1 | Correct', 'Color', color_type{j});
                    elseif j == 4
                        ylabel('Blocktype 1 | Error', 'Color', color_type{j});
                    end
                end
                ph = patch([t_before_tone_onset t_before_tone_onset+300 t_before_tone_onset+300 t_before_tone_onset],[0 0 ntrials+0.5 ntrials+0.5],...
                [.8 .8 .8],'edgecolor', [.8 .8 .8]);
                alpha(ph,.4);                              
            end
        end
        
        j = 5;
        subplot(5,length(fre_type), ifre + (j-1) * length(fre_type));
        hold on
        x = [1:numbin]*bin_length;
        for j5 =1:4
            if sum(fre_corr_block_idx{ifre,j5}) >= 3
                y(:,j5) = smooth(mean(squeeze(np2tone(iunit, fre_corr_block_idx{ifre,j5}, :)),1, 'omitnan'), 40);
                plot(x, y(:,j5), 'Color', color_type{j5});
            end
        end
        
        ymax = max(max(y(10:end-10,:)));
        ph2 = patch([t_before_tone_onset t_before_tone_onset+300 t_before_tone_onset+300 t_before_tone_onset],[0 0 ymax+1 ymax+1],...
        [.8 .8 .8],'edgecolor',[.8 .8 .8]);
    
        alpha(ph2,.4);
        
        xlim([0,bin_length*size(np2tone,3)]);
        ylim([0, ymax+1]);
        y = [];
        if ifre == 1
            ylabel('Firing rate (Hz)');
            xlabel('Time (ms)');
        end
        hold off;
        
    end
%     saveas(gcf, [fullfile(save_rasterplot_folder, savename) '.png']);
%     saveas(gcf, [fullfile(save_rasterplot_folder, savename) '.fig']);
%    fprintf([num2str(iunit), '/',num2str(size(np_data_cell, 2)), 'saved.\n']);
%    close;
end
