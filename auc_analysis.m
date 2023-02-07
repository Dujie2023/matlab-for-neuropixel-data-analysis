%% cal auc for block, corr/error, left/right choice
% np_data_pathway = 'G:\NPrecording_XY\Cat_A20220101_b107a07_NPSess06_g0_imec1'; % %



neuroactivity_bftone = [];  neuroactivity_aftone = [];
neuroactivity_bfanswer = [];    neuroactivity_afanswer = [];
for icluster = 1:size(np2tone,1)
    ineuroactivity_bftone = mean(np2tone(icluster,:,1:bin_stim_onset-1), 3);
    [auc_block_bftone(icluster),auc_block_bftone1(icluster)] = cal_auc2(ineuroactivity_bftone, blocktype_use);
    [auc_corr_bftone(icluster),auc_corr_bftone1(icluster)] = cal_auc2(ineuroactivity_bftone, TriCorr_use);
    [auc_left_bftone(icluster),auc_left_bftone1(icluster)] = cal_auc2(ineuroactivity_bftone, TriLeft_use);
    neuroactivity_bftone = [neuroactivity_bftone; ineuroactivity_bftone];
    
    ineuroactivity_aftone = mean(np2tone(icluster,:,bin_stim_onset:200), 3);    
    [auc_block_aftone(icluster),auc_block_aftone1(icluster)] = cal_auc2(ineuroactivity_aftone, blocktype_use);
    [auc_corr_aftone(icluster),auc_corr_aftone1(icluster)] = cal_auc2(ineuroactivity_aftone, TriCorr_use);
    [auc_left_aftone(icluster),auc_left_aftone1(icluster)] = cal_auc2(ineuroactivity_aftone, TriLeft_use);
    neuroactivity_aftone = [neuroactivity_aftone; ineuroactivity_aftone];
    
    ineuroactivity_bfanswer = mean(np2answer(icluster,:,1:bin_answer-1), 3);
    [auc_block_bfanswer(icluster),auc_block_bfanswer1(icluster)] = cal_auc2(ineuroactivity_bfanswer, blocktype_use);
    [auc_corr_bfanswer(icluster),auc_corr_bfanswer1(icluster)] = cal_auc2(ineuroactivity_bfanswer, TriCorr_use);
    [auc_left_bfanswer(icluster),auc_left_bfanswer1(icluster)] = cal_auc2(ineuroactivity_bfanswer, TriLeft_use);
    neuroactivity_bfanswer = [neuroactivity_bfanswer; ineuroactivity_bfanswer];
    
    ineuroactivity_afanswer = mean(np2answer(icluster,:,bin_answer:end), 3);    
    [auc_block_afanswer(icluster),auc_block_afanswer1(icluster)] = cal_auc2(ineuroactivity_afanswer, blocktype_use);
    [auc_corr_afanswer(icluster),auc_corr_afanswer1(icluster)] = cal_auc2(ineuroactivity_afanswer, TriCorr_use);
    [auc_left_afanswer(icluster),auc_left_afanswer1(icluster)] = cal_auc2(ineuroactivity_afanswer, TriLeft_use);
    neuroactivity_afanswer = [neuroactivity_afanswer; ineuroactivity_afanswer];

end

SI_block = 2*([auc_block_bftone1; auc_block_aftone1; auc_block_bfanswer1;auc_block_afanswer1]-0.5);
SI_corr = 2*([auc_corr_bftone1; auc_corr_aftone1; auc_corr_bfanswer1;auc_corr_afanswer1]-0.5);
SI_left = 2*([auc_left_bftone1; auc_left_aftone1; auc_left_bfanswer1;auc_left_afanswer1]-0.5);
AUC_block = [auc_block_bftone; auc_block_aftone; auc_block_bfanswer;auc_block_afanswer];
AUC_corr = [auc_corr_bftone; auc_corr_aftone; auc_corr_bfanswer;auc_corr_afanswer];
AUC_left = [auc_left_bftone; auc_left_aftone; auc_left_bfanswer;auc_left_afanswer];

save(fullfile(np_data_pathway, 'auc'), 'AUC_block','AUC_corr','AUC_left','neuroactivity_bftone',...
    'neuroactivity_aftone','neuroactivity_bfanswer', 'neuroactivity_afanswer','blocktype_use');
%% cal auc significance threshold
all_neuroactivity = neuroactivity_bftone;
[block_sig_thre05, block_sig_thre01] = auc_sigthreshold(all_neuroactivity, blocktype_use, 1000);
[corr_sig_thre05, corr_sig_thre01] = auc_sigthreshold(all_neuroactivity, TriCorr_use, 1000);
[left_sig_thre05, left_sig_thre01] = auc_sigthreshold(all_neuroactivity, TriLeft_use, 1000);
save(fullfile(np_data_pathway, 'sig_threshold'), 'block_sig_thre05','block_sig_thre01','corr_sig_thre05','corr_sig_thre01',...
    'left_sig_thre05','left_sig_thre01');
load(fullfile(np_data_pathway, 'sig_threshold'));

%% plot for block
SigFracBlock_areas = nan(length(Abs),4);%fraction of neurons signicantly discriminate block in different areas

%aline 2 answer
figure;
hold on
for iAb = 1:length(Abs)
    SigFracBlock_areas(iAb,3) = sum(AUC_block(3,Ab_idx{iAb})>block_sig_thre01)/length(Ab_idx{iAb});
    SigFracBlock_areas(iAb,4) = sum(AUC_block(4,Ab_idx{iAb})>block_sig_thre01)/length(Ab_idx{iAb});
    scatter(SI_block(3,Ab_idx{iAb}),SI_block(4,Ab_idx{iAb}),20,'filled');    
end
limset = max(max(xlim, ylim));
xlim([-1*limset limset]);
ylim([-1*limset limset]);
xlabel('Selectivity index before answer');
ylabel('Selectivity after answer');
title('Selectivity Index for Block');

xdiag = [-1*limset:0.05:limset];
ydiag = xdiag;
plot(xdiag,ydiag,'k--');

legend(usedAb,'Location','Best');
legend('boxoff');

[r, p] = corr(SI_block(3,:)',SI_block(4,:)');
txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
text(360,320, txtr, 'Units', 'pixels');
text(360,300, txtp, 'Units', 'pixels');
hold off;
saveas(gcf, [np_data_pathway, '\auc\scatter_block_answer_area.png']);
saveas(gcf, [np_data_pathway, '\auc\scatter_block_answer_area.fig']);
close;

%aline2tone
figure;
hold on
for iAb = 1:length(Abs)
    SigFracBlock_areas(iAb,1) = sum(AUC_block(1,Ab_idx{iAb})>block_sig_thre01)/length(Ab_idx{iAb});
    SigFracBlock_areas(iAb,2) = sum(AUC_block(2,Ab_idx{iAb})>block_sig_thre01)/length(Ab_idx{iAb});
    scatter(SI_block(1,Ab_idx{iAb}),SI_block(2,Ab_idx{iAb}),20,'filled');    
end
limset = max(max(xlim, ylim));
xlim([-1*limset limset]);
ylim([-1*limset limset]);
xlabel('Selectivity index before tone');
ylabel('Selectivity after tone');
title('Selectivity Index for Block');

xdiag = [-1*limset:0.05:limset];
ydiag = xdiag;
plot(xdiag,ydiag,'k--');

legend(usedAb,'Location','Best');
legend('boxoff');

[r, p] = corr(SI_block(1,:)',SI_block(2,:)');
txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
text(360,320, txtr, 'Units', 'pixels');
text(360,300, txtp, 'Units', 'pixels');
hold off;
saveas(gcf, [np_data_pathway, '\auc\scatter_block_tone_area.png']);
saveas(gcf, [np_data_pathway, '\auc\scatter_block_tone_area.fig']);
close;

%show fraction of neurons can discriminate
figure;
%set(gcf, 'Position', [100 100 700 600]);
hold on
x = 1:4;
y = SigFracBlock_areas;
for i =1:length(Abs)
    plot(x,y(i,:),'.-','LineWidth',2, 'markersize',20);
end
legend(usedAb,'Location','best','FontSize',12);
legend('boxoff');

ylabel('Fraction of units','FontSize',20);
title('Fraction of units significantly discriminate block');
set(gca,'xTick',[1:4],'xTickLabel',{'bftone','aftone','bfanswer','afanswer'}, 'FontSize',10);

saveas(gcf,fullfile(np_data_pathway, 'auc',  'block_perSI.png'));
saveas(gcf,fullfile(np_data_pathway, 'auc',  'block_perSI.fig'));
close;

%colorplot
auc1block2tone = [];
auc1block2answer = [];
ave_nbin =10;
for iunit = 1:size(np2answer,1) 
    for jbin = 1:size(np2tone,3)/ave_nbin
        temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
        [~, auc1block2tone(iunit,jbin)]=cal_auc2(temp, blocktype_use);     
    end
    for jbin = 1:size(np2answer,3)/ave_nbin
        temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
        [~, auc1block2answer(iunit,jbin)]=cal_auc2(temp, blocktype_use);        
    end
end

SIblock2tone = 2*(auc1block2tone-0.5);
SIblock2answer = 2*(auc1block2answer-0.5);

SIblock2tone_sorted = [];
SIblock2answer_sorted = [];
SIblock2tone(isnan(SIblock2tone(:))) = 0;
SIblock2answer(isnan(SIblock2answer(:))) = 0;
SIblock_mean = mean(SIblock2answer,2);
[~, sort_idx] = sort(SIblock_mean);
for i = 1:length(SIblock_mean)
    SIblock2answer_sorted(i,:) = SIblock2answer(sort_idx(i),:);
    SIblock2tone_sorted(i,:) = SIblock2tone(sort_idx(i),:);
end

figure; 
set(gcf, 'Position', [100 100 1400 700]);
subplot(1,2,2);
imagesc(SIblock2answer_sorted);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
colorbar;
temp1 = caxis;
set(gca, 'Xtick', [0:10:size(np2answer,3)/ave_nbin], 'XTicklabel', ...
    100*[0:10:size(np2answer,3)/ave_nbin]);
title(['Selectivity Index for Block ',' 2answer'], 'FontSize', 14);
%answer time
hold on;
xanswer = bin_answer*ones(1,size(np2answer,1))/ave_nbin;
yanswer = [1:size(np2answer,1)];
line(xanswer, yanswer, 'LineWidth', 1.5,'LineStyle', '--', 'Color', 'k');
hold off;

subplot(1,2,1) 
imagesc(SIblock2tone_sorted(:,1:2000/(bin_length*ave_nbin)));
caxis(temp1);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
set(gca, 'Xtick', [0:10:2000/(bin_length*ave_nbin)], 'XTicklabel', ...
    [0:10*(bin_length*ave_nbin):2000]);
title(['Selectivity Index for Block ',' 2tone'],'FontSize', 14);
%tone time
ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin, ...
      [size(SIblock2tone, 1) 0 0  size(SIblock2tone, 1)], 'w', 'edgecolor', 'none');
alpha(ph,.4);
filename = fullfile(np_data_pathway, 'auc', 'block_allunits');
saveas(gcf,[filename  '.png']);
saveas(gcf,[filename  '.fig']);

%% plot for corr/error
SigFracCorr_areas = nan(length(Abs),4);%fraction of neurons can discriminate corr in different areas

figure;
hold on
for iAb = 1:length(Abs)
    SigFracCorr_areas(iAb,3) = sum(AUC_corr(3,Ab_idx{iAb})>corr_sig_thre01)/length(Ab_idx{iAb});
    SigFracCorr_areas(iAb,4) = sum(AUC_corr(4,Ab_idx{iAb})>corr_sig_thre01)/length(Ab_idx{iAb});
    scatter(SI_corr(3,Ab_idx{iAb}),SI_corr(4,Ab_idx{iAb}),20,'filled');    
end
limset = max(max(xlim, ylim));
xlim([-1*limset limset]);
ylim([-1*limset limset]);
xlabel('Selectivity index before answer');
ylabel('Selectivity after answer');
title('Selectivity Index for Corr/Error');

xdiag = [-1*limset:0.05:limset];
ydiag = xdiag;
plot(xdiag,ydiag,'k--');

legend(usedAb,'Location','Best');
legend('boxoff');

[r, p] = corr(SI_corr(3,:)',SI_corr(4,:)');
txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
text(360,320, txtr, 'Units', 'pixels');
text(360,300, txtp, 'Units', 'pixels');
hold off;
saveas(gcf, [np_data_pathway, '\auc\scatter_corr_answer_area.png']);
saveas(gcf, [np_data_pathway, '\auc\scatter_corr_answer_area.fig']);
close;

%aline2tone
figure;
hold on
for iAb = 1:length(Abs)
    SigFracCorr_areas(iAb,1) = sum(AUC_corr(1,Ab_idx{iAb})>corr_sig_thre01)/length(Ab_idx{iAb});
    SigFracCorr_areas(iAb,2) = sum(AUC_corr(2,Ab_idx{iAb})>corr_sig_thre01)/length(Ab_idx{iAb});
    scatter(SI_corr(1,Ab_idx{iAb}),SI_corr(2,Ab_idx{iAb}),20,'filled');    
end
limset = max(max(xlim, ylim));
xlim([-1*limset limset]);
ylim([-1*limset limset]);
xlabel('Selectivity index before tone');
ylabel('Selectivity after tone');
title('Selectivity Index for Corr/Erorr');

xdiag = [-1*limset:0.05:limset];
ydiag = xdiag;
plot(xdiag,ydiag,'k--');

legend(usedAb,'Location','Best');
legend('boxoff');

[r, p] = corr(SI_corr(1,:)',SI_corr(2,:)');
txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
text(360,320, txtr, 'Units', 'pixels');
text(360,300, txtp, 'Units', 'pixels');
hold off;
saveas(gcf, [np_data_pathway, '\auc\scatter_corr_tone_area.png']);
saveas(gcf, [np_data_pathway, '\auc\scatter_corr_tone_area.fig']);
close;

%show fraction of neurons can discriminate
figure;
%set(gcf, 'Position', [100 100 700 600]);
hold on
x = 1:4;
y = SigFracCorr_areas;
for i =1:length(Abs)
    plot(x,y(i,:),'.-','LineWidth',2, 'markersize',20);
end
legend(usedAb,'Location','best','FontSize',12);
legend('boxoff');

ylabel('Fraction of units','FontSize',20);
title('Fraction of units significantly discriminate corr/error');
set(gca,'xTick',[1:4],'xTickLabel',{'bftone','aftone','bfanswer','afanswer'}, 'FontSize',10);

saveas(gcf,fullfile(np_data_pathway, 'auc',  'corr_perSI.png'));
saveas(gcf,fullfile(np_data_pathway, 'auc',  'corr_perSI.fig'));
close;

%color plot
ave_nbin =10;
auc1corr2tone = [];
auc1corr2answer = [];
for iunit = 1:size(np2answer,1) 
    for jbin = 1:size(np2tone,3)/ave_nbin
        temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
        [~, auc1corr2tone(iunit,jbin)]=cal_auc2(temp, TriCorr_use);     
    end
    for jbin = 1:size(np2answer,3)/ave_nbin
        temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
        [~, auc1corr2answer(iunit,jbin)]=cal_auc2(temp, TriCorr_use);        
    end
end

SIcorr2tone = 2*(auc1corr2tone-0.5);
SIcorr2answer = 2*(auc1corr2answer-0.5);

SIcorr2tone_sorted = [];
SIcorr2answer_sorted = [];
SIcorr2tone(isnan(SIcorr2tone(:))) = 0;
SIcorr2answer(isnan(SIcorr2answer(:))) = 0;
SIcorr_mean = mean(SIcorr2answer,2);
[~, sort_idx] = sort(SIcorr_mean);
for i = 1:length(SIcorr_mean)
    SIcorr2answer_sorted(i,:) = SIcorr2answer(sort_idx(i),:);
    SIcorr2tone_sorted(i,:) = SIcorr2tone(sort_idx(i),:);
end

figure; 
set(gcf, 'Position', [100 100 1400 700]);
subplot(1,2,2);
imagesc(SIcorr2answer_sorted);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
colorbar;
temp1 = caxis;
set(gca, 'Xtick', [0:10:size(np2answer,3)/ave_nbin], 'XTicklabel', ...
    100*[0:10:size(np2answer,3)/ave_nbin]);
title(['Selectivity Index for Corr/Error ',' 2answer'], 'FontSize', 14);
%answer time
hold on;
xanswer = bin_answer*ones(1,size(np2answer,1))/ave_nbin;
yanswer = [1:size(np2answer,1)];
line(xanswer, yanswer, 'LineWidth', 1.5,'LineStyle', '--', 'Color', 'k');
hold off;

subplot(1,2,1) 
imagesc(SIcorr2tone_sorted(:,1:2000/(bin_length*ave_nbin)));
caxis(temp1);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
set(gca, 'Xtick', [0:10:2000/(bin_length*ave_nbin)], 'XTicklabel', ...
    [0:10*(bin_length*ave_nbin):2000]);
title(['Selectivity Index for Corr/Error ',' 2tone'],'FontSize', 14);
%tone time
ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin, ...
      [size(SIcorr2tone, 1) 0 0  size(SIcorr2tone, 1)], 'w', 'edgecolor', 'none');
alpha(ph,.4);
filename = fullfile(np_data_pathway, 'auc', 'corr_allunits');
saveas(gcf,[filename  '.png']);
saveas(gcf,[filename  '.fig']);
%% plot for left/right

SigFracLeft_areas = nan(length(Abs),4);%fraction of neurons can discriminate in different areas

figure;
hold on
for iAb = 1:length(Abs)
    SigFracLeft_areas(iAb,3) = sum(AUC_left(3,Ab_idx{iAb})>left_sig_thre01)/length(Ab_idx{iAb});
    SigFracLeft_areas(iAb,4) = sum(AUC_left(4,Ab_idx{iAb})>left_sig_thre01)/length(Ab_idx{iAb});
    scatter(SI_left(3,Ab_idx{iAb}),SI_left(4,Ab_idx{iAb}),20,'filled');    
end
limset = max(max(xlim, ylim));
xlim([-1*limset limset]);
ylim([-1*limset limset]);
xlabel('Selectivity index before answer');
ylabel('Selectivity after answer');
title('Selectivity Index for left choice');

xdiag = [-1*limset:0.05:limset];
ydiag = xdiag;
plot(xdiag,ydiag,'k--');

legend(usedAb,'Location','Best');
legend('boxoff');

[r, p] = corr(SI_left(3,:)',SI_left(4,:)');
txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
text(360,320, txtr, 'Units', 'pixels');
text(360,300, txtp, 'Units', 'pixels');
hold off;
saveas(gcf, [np_data_pathway, '\auc\scatter_left_answer_area.png']);
saveas(gcf, [np_data_pathway, '\auc\scatter_left_answer_area.fig']);
close;

%aline2tone
figure;
hold on
for iAb = 1:length(Abs)
    SigFracLeft_areas(iAb,1) = sum(AUC_left(1,Ab_idx{iAb})>left_sig_thre01)/length(Ab_idx{iAb});
    SigFracLeft_areas(iAb,2) = sum(AUC_left(2,Ab_idx{iAb})>left_sig_thre01)/length(Ab_idx{iAb});
    scatter(SI_left(1,Ab_idx{iAb}),SI_left(2,Ab_idx{iAb}),20,'filled');    
end
limset = max(max(xlim, ylim));
xlim([-1*limset limset]);
ylim([-1*limset limset]);
xlabel('Selectivity index before tone');
ylabel('Selectivity after tone');
title('Selectivity Index for left choice');

xdiag = [-1*limset:0.05:limset];
ydiag = xdiag;
plot(xdiag,ydiag,'k--');

legend(usedAb,'Location','Best');
legend('boxoff');

[r, p] = corr(SI_left(1,:)',SI_left(2,:)');
txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
text(360,320, txtr, 'Units', 'pixels');
text(360,300, txtp, 'Units', 'pixels');
hold off;
saveas(gcf, [np_data_pathway, '\auc\scatter_left_tone_area.png']);
saveas(gcf, [np_data_pathway, '\auc\scatter_left_tone_area.fig']);
close;

%show fraction of neurons can discriminate
figure;
%set(gcf, 'Position', [100 100 700 600]);
hold on
x = 1:4;
y = SigFracLeft_areas;
for i =1:length(Abs)
    plot(x,y(i,:),'.-','LineWidth',2, 'markersize',20);
end
legend(usedAb,'Location','best','FontSize',12);
legend('boxoff');

ylabel('Fraction of units','FontSize',20);
title('Fraction of units significantly discriminate left choice');
set(gca,'xTick',[1:4],'xTickLabel',{'bftone','aftone','bfanswer','afanswer'}, 'FontSize',10);

saveas(gcf,fullfile(np_data_pathway, 'auc',  'left_perSI.png'));
saveas(gcf,fullfile(np_data_pathway, 'auc',  'left_perSI.fig'));

%color plot
ave_nbin =10;
auc1left2tone = [];
auc1left2answer = [];
for iunit = 1:size(np2answer,1) 
    for jbin = 1:size(np2tone,3)/ave_nbin
        temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
        [~, auc1left2tone(iunit,jbin)]=cal_auc2(temp, TriLeft_use);     
    end
    for jbin = 1:size(np2answer,3)/ave_nbin
        temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
        [~, auc1left2answer(iunit,jbin)]=cal_auc2(temp, TriLeft_use);        
    end
end

SIleft2tone = 2*(auc1left2tone-0.5);
SIleft2answer = 2*(auc1left2answer-0.5);

SIleft2tone_sorted = [];
SIleft2answer_sorted = [];
SIleft2tone(isnan(SIleft2tone(:))) = 0;
SIleft2answer(isnan(SIleft2answer(:))) = 0;
SIleft_mean = mean(SIleft2answer,2);
[~, sort_idx] = sort(SIleft_mean);
for i = 1:length(SIleft_mean)
    SIleft2answer_sorted(i,:) = SIleft2answer(sort_idx(i),:);
    SIleft2tone_sorted(i,:) = SIleft2tone(sort_idx(i),:);
end

figure; 
set(gcf, 'Position', [100 100 1400 700]);
subplot(1,2,2);
imagesc(SIleft2answer_sorted);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
colorbar;
temp1 = caxis;
set(gca, 'Xtick', [0:10:size(np2answer,3)/ave_nbin], 'XTicklabel', ...
    100*[0:10:size(np2answer,3)/ave_nbin]);
title(['Selectivity Index for Corr/Error ',' 2answer'], 'FontSize', 14);
%answer time
hold on;
xanswer = bin_answer*ones(1,size(np2answer,1))/ave_nbin;
yanswer = [1:size(np2answer,1)];
line(xanswer, yanswer, 'LineWidth', 1.5,'LineStyle', '--', 'Color', 'k');
hold off;

subplot(1,2,1) 
imagesc(SIleft2tone_sorted(:,1:2000/(bin_length*ave_nbin)));
caxis(temp1);
ylabel('Unit Number','FontSize',12);
xlabel('Time (ms)','FontSize',12);
set(gca, 'Xtick', [0:10:2000/(bin_length*ave_nbin)], 'XTicklabel', ...
    [0:10*(bin_length*ave_nbin):2000]);
title(['Selectivity Index for Corr/Error ',' 2tone'],'FontSize', 14);
%tone time
ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin, ...
      [size(SIleft2tone, 1) 0 0  size(SIleft2tone, 1)], 'w', 'edgecolor', 'none');
alpha(ph,.4);
filename = fullfile(np_data_pathway, 'auc', 'left_allunits');
saveas(gcf,[filename  '.png']);
saveas(gcf,[filename  '.fig']);

%% plot summary for different brain areas
numSig01Areas = nan(3,length(Abs),4);
varSI = {'Block','Corr/Error','Left/Right'};
SIthreshold = 0.5;
for iAb = 1:length(Abs)
     numSig01Areas(1, iAb, :) = sum(AUC_block(:,Ab_idx{iAb})>block_sig_thre01, 2);
     numSig01Areas(2, iAb, :) = sum(AUC_corr(:,Ab_idx{iAb})>corr_sig_thre01, 2);
     numSig01Areas(3, iAb, :) = sum(AUC_left(:,Ab_idx{iAb})>left_sig_thre01, 2);
end
for iAb = 1:length(Abs)
    
    figure;
    hold on
    x = 1:4;
    y = squeeze(numSig01Areas(:, iAb,:));
    for i =1:length(varSI)
        plot(x,y(i,:),'.-','LineWidth',2, 'markersize',20);
    end
    legend(varSI,'Location','best','FontSize',12);
    legend('boxoff');

    ylabel('Number of units','FontSize',20);
    title(usedAb{iAb});
    set(gca,'xTick',[1:4],'xTickLabel',{'bftone','aftone','bfanswer','afanswer'}, 'FontSize',10);

    saveas(gcf,fullfile(np_data_pathway, 'auc',  [usedAb{iAb} '_Sig01.png']));
    saveas(gcf,fullfile(np_data_pathway, 'auc',  [usedAb{iAb} '_Sig01.fig']));
end
% %% plot venn
% % 3sets: block/corr/left
% % 4sub plots for different time windows
% % plot for different brain area
% for iAb = 5:6%1:length(Abs)
%     
%     for timewindow = 1:4
%         tempblock{timewindow,iAb} = find(abs(SI_block(timewindow,Ab_idx{iAb}))>0.5);    
%         tempcorr{timewindow,iAb} = find(abs(SI_corr(timewindow,Ab_idx{iAb}))>0.5);
%         templeft{timewindow,iAb} = find(abs(SI_left(timewindow,Ab_idx{iAb}))>0.5);
%         
%     end
%     
%     A = [length(tempblock{timewindow,iAb}), length(tempcorr{timewindow,iAb}), length(templeft{timewindow,iAb})];
%     I = [length(intersect(tempblock{timewindow,iAb},tempcorr{timewindow,iAb})), length(intersect(tempblock{timewindow,iAb},templeft{timewindow,iAb})),...
%          length(intersect(tempcorr{timewindow,iAb},templeft{timewindow,iAb})),length(intersect(intersect(tempblock{timewindow,iAb}, tempcorr{timewindow,iAb}),templeft{timewindow,iAb}))];
% 
%     figure
%     venn(A,I,'Plot','off');
% end