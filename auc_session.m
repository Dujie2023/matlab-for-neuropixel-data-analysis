function auc_session(np_data_pathway, cur_label,labelname,nbin_aftone, nbin_afanswer,useUtest)    
%% data loading
%cur_label: label vector(ntrials), e.g Choice, Outcome, Block
%label name: char
%------time window setting
%nbin_aftone: how many bins after stime tone onset are used, 10ms for a bin;
%nbin_afanswer: how many bins after answer are used, 10ms for a bin;
%useUtest: choose algorithm for calculating, useUtest=1: U test
save_folder = fullfile(np_data_pathway, 'auc');
if ~exist(save_folder)
    mkdir(save_folder);
end

load(fullfile(np_data_pathway, 'np_extracted.mat'));
%if unit number is 0, skip analysis
if isempty(Spike_firingrate_bin)
    warning('The good unit number of this probe is 0.');
    return;
end

np2tone = Spike_firingrate_bin;%neural activity, alined to tone onset time, matrix(nunits*ntrials*nbins)
np2answer =Spike_firingrate_bin_outcome; %neural activity, alined to answer time, matrix(nunits*ntrials*nbins)
neuroactivity_bftone = mean(np2tone(:,:,1:bin_stim_onset-1), 3);
neuroactivity_aftone = mean(np2tone(:,:,bin_stim_onset:bin_stim_onset+nbin_aftone), 3);    
neuroactivity_bfanswer = mean(np2answer(:,:,1:bin_answer-1), 3);
neuroactivity_afanswer = mean(np2answer(:,:,bin_answer:bin_answer+nbin_afanswer), 3);    
%% cal auc 
auc_bftone = nan(size(neuroactivity_bftone,1),1);
auc_aftone = nan(size(neuroactivity_bftone,1),1);
auc_bfanswer = nan(size(neuroactivity_bftone,1),1);
auc_afanswer = nan(size(neuroactivity_bftone,1),1);

if useUtest == 1
    for icluster = 1:size(neuroactivity_bftone,1)
        auc_bftone(icluster) = cal_auc_Utest(neuroactivity_bftone(icluster,:), cur_label);
        auc_aftone(icluster) = cal_auc_Utest(neuroactivity_aftone(icluster,:), cur_label);
        auc_bfanswer(icluster) = cal_auc_Utest(neuroactivity_bfanswer(icluster,:), cur_label);
        auc_afanswer(icluster) = cal_auc_Utest(neuroactivity_afanswer(icluster,:), cur_label);
    end
else
   for icluster = 1:size(neuroactivity_bftone,1)
        auc_bftone(icluster) = cal_auc2(neuroactivity_bftone(icluster,:), cur_label);
        auc_aftone(icluster) = cal_auc2(neuroactivity_aftone(icluster,:), cur_label);
        auc_bfanswer(icluster) = cal_auc2(neuroactivity_bfanswer(icluster,:), cur_label);
        auc_afanswer(icluster) = cal_auc2(neuroactivity_afanswer(icluster,:), cur_label);
    end
end

% %plot for the session, neurons from different areas will not be
% %discriminated
figure;
scatter(auc_bftone,  auc_aftone);
xlabel(['Before Tone AUC ',labelname]);
ylabel(['After Tone AUC ',labelname]);
saveas(gcf, fullfile(save_folder,[labelname 'bftone_vs_aftone_auc.png']));
saveas(gcf, fullfile(save_folder,[labelname 'bftone_vs_aftone_auc.fig']));
close;

figure;
scatter(auc_bfanswer,  auc_afanswer);
xlabel(['Before answer AUC ',labelname]);
ylabel(['After answer AUC ',labelname]);
saveas(gcf, fullfile(save_folder,[labelname, 'bfanswer_vs_afafanswer_auc.png']));
saveas(gcf, fullfile(save_folder,[labelname, 'bfanswer_vs_afafanswer_auc.fig']));
close;

%% cal threshold of significance
[sig_thre05_aline2tone, sig_thre01_aline2tone] = parfor_auc_sigthreshold([neuroactivity_bftone;neuroactivity_aftone], cur_label, 1000, save_folder, labelname, useUtest);
[sig_thre05_aline2answer, sig_thre01_aline2answer] = parfor_auc_sigthreshold([neuroactivity_bfanswer;neuroactivity_afanswer], cur_label, 1000, save_folder, labelname, useUtest);

%% save results 
save(fullfile(save_folder,[labelname,' auc.mat']), 'auc_bftone', 'auc_aftone',...
  'auc_bfanswer', 'auc_afanswer','sig_thre05_aline2tone', 'sig_thre01_aline2tone',...
  'sig_thre05_aline2answer','sig_thre01_aline2answer','cur_label','labelname'...
  ,'nbin_aftone', 'nbin_afanswer');

%% plot for small time window
% ave_nbin_auc = 20;
% numbin2tone = floor(size(np2tone,3)/ave_nbin_auc);
% numbin2answer = floor(size(np2answer,3)/ave_nbin_auc);
% num_unit = size(np2answer,1);
% auc2tone = nan(num_unit,numbin2tone);
% auc2answer = nan(num_unit,numbin2answer);
% 
% for iunit = 1:size(np2answer,1) 
%     for jbin = 1:numbin2tone
%         temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin_auc-ave_nbin_auc+1:jbin*ave_nbin_auc), 3));
%         auc2tone(iunit,jbin) = cal_auc2(temp, cur_label);        
%     end
%     for jbin = 1:numbin2answer
%         temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin_auc-ave_nbin_auc+1:jbin*ave_nbin_auc), 3));
%         auc2answer(iunit,jbin) = cal_auc2(temp, cur_label);        
%     end
% end
% 
% %sort auc matrix
% auc_mean = mean(auc2answer,2);
% [~, sort_idx] = sort(auc_mean);
% auc2answer_sorted = auc2answer(sort_idx,:);
% auc2tone_sorted = auc2tone(sort_idx,:);
% 
% figure; 
% set(gcf, 'Position', [100 100 1400 700]);
% subplot(1,2,2);
% imagesc(auc2answer_sorted);
% ylabel('Unit Number','FontSize',12);
% xlabel('Time (ms)','FontSize',12);
% colorbar;
% temp1 = caxis;
% set(gca, 'Xtick', [0:1000/(bin_length*ave_nbin_auc):size(np2answer,3)/ave_nbin_auc], 'XTicklabel', ...
%    [0:1000:size(np2answer,3)*bin_length]);
% title([labelname ' AUC ', ' 2answer'], 'FontSize', 14);
% %answer time
% hold on;
% xanswer = bin_answer*ones(1,size(np2answer,1))/ave_nbin_auc;
% yanswer = [1:size(np2answer,1)];
% line(xanswer, yanswer, 'LineWidth', 1.5,'LineStyle', '--', 'Color', 'k');
% hold off;
% 
% subplot(1,2,1) 
% imagesc(auc2tone_sorted(:,1:2000/(bin_length*ave_nbin_auc)));
% caxis(temp1);
% ylabel('Unit Number','FontSize',12);
% xlabel('Time (ms)','FontSize',12);
% set(gca, 'Xtick', [0:1000/(bin_length*ave_nbin_auc):2000/(bin_length*ave_nbin_auc)], 'XTicklabel', ...
%     [0:1000:2000]);
% title([labelname ' AUC ',' 2tone'],'FontSize', 14);
% %tone time
% ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin_auc, ...
%       [size(auc2tone, 1) 0 0  size(auc2tone, 1)], 'w', 'edgecolor', 'none');
% alpha(ph,.4);
% filename = fullfile(save_folder,[labelname '_allunit']);
% saveas(gcf,[filename  '.png']);
% saveas(gcf,[filename  '.fig']);

%% plot example neuron (only for block label)
if strcmp(labelname, 'Block')

    block = double(cur_label);
    blockswitchtri = find(diff(block));
    [auc_block_bftone_sorted,example_idx] = sort(auc_bftone,'descend');

    load(fullfile(np_data_pathway,'RL_model','BSRL.mat'));
    imodelvar = P_bound_low(1:end-1);

    load(fullfile(np_data_pathway, 'probeNPSess.mat'));
    ChannelAreaStrs = probeNPSess.ChannelAreaStrs{1,2};

    ienum = 20;%the number of example neurons
    if size(np2tone,1)<20
        ienum = size(np2tone,1);
    end

    for ie = 1:ienum
        example_num = example_idx(ie);
        iechannel = probeNPSess.ChannelUseds_id(example_num);
        ieAreaStr = ChannelAreaStrs{iechannel};
        figure;
        x = 1:length(imodelvar);
        y1 = smooth(imodelvar,10);
        y2 = smooth(mean(squeeze(np2answer(example_num,:,1:(bin_stim_onset-1))),2),10);
        plot(x,y1,'k','LineWidth',1);

        hold on;
        axisrange = axis;
        yrange = axisrange(3:4);
        for iblockswitch = 1:length(blockswitchtri)
            plot(ones(1,20)*blockswitchtri(iblockswitch),linspace(yrange(1),yrange(2),20),'k--','LineWidth',1');
        end

        yyaxis right
        plot(x,y2,'b','LineWidth',1);
        set(gca,'ycolor','b');
        ylabel('Firing rate (Hz)');
        set(gca,'TickDir','out');

        yyaxis left
        ylabel( 'P boundary low(t)','FontSize',12);
        xlabel('# Trial','FontSize',12);
        set(gca,'TickDir','out');

        auc_ie = auc_block_bftone_sorted(ie);
        title(['unit ' num2str(example_num) ' (Area: ' ieAreaStr ')' ' AUC = ' num2str(auc_ie,2)]);

        box off;
        hold off;

        saveas(gcf,[save_folder, '\example unit ', num2str(example_num),'.png']);
        close;

    end

end


end