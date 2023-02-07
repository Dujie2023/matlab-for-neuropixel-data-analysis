for iunit = 1:size(np_data_cell, 2)
    
    figure;
    subplot(5,7,1);

    
    for t = 1:ntrials
    t_unit_spike = np_data_cell{t,iunit}; % Spike timings in the t_th trial for unit_idx
    nspikes = numel(t_unit_spike); % number of spikes
        for ii = 1:nspikes % for every spike
            pl = line([t_unit_spike(ii) t_unit_spike(ii)],[t-0.5 t+0.5] );
            pl.Color = 'Black';% draw a black vertical line of length 1 at time t(x) and at trial t (y)
        end
    end
    xlabel('Time (ms)'); % Time is in milliseconds
    ylabel('Trial number ');
    ylim([0,ntrials]);
    xlim([0,8000]);
    ph = patch([t_before_tone_onset t_before_tone_onset+300 t_before_tone_onset+300 t_before_tone_onset],[0 0 ntrials+0.5 ntrials+0.5],...
        [.8 .8 .8],'edgecolor', [.8 .8 .8]);
    alpha(ph,.3);
    savename = ['Cluster ',num2str(Used_ID(iunit))];
    title(savename);

    subplot(2,1,2);
    x = [1:numbin]*bin_length;
    y = smooth(mean(squeeze(np_data(iunit, :, :)),1, 'omitnan')*1000, 5);
    
    plot(x,y,'k');
    ph2 = patch([t_before_tone_onset t_before_tone_onset+300 t_before_tone_onset+300 t_before_tone_onset],[0 0 max(y)+1 max(y)+1],...
        [.8 .8 .8],'edgecolor',[.8 .8 .8]);
    alpha(ph2,.3);
    
    xlim([0,8000]);
    ylim([0,max(y)+1]);
    ylabel('Firing rate (Hz)');
    xlabel('Time (ms)');

    saveas(gcf,[savename  '.png']);