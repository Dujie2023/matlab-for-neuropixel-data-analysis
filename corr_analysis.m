%% cal all R,P
R_allvar_bftone = zeros(size(modelvar,1),size(np2tone,1));
R_allvar_aftone = zeros(size(modelvar,1),size(np2tone,1));
R_allvar_bfanswer = zeros(size(modelvar,1),size(np2tone,1));
R_allvar_afanswer = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_bftone = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_aftone = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_bfanswer = zeros(size(modelvar,1),size(np2tone,1));
P_allvar_afanswer = zeros(size(modelvar,1),size(np2tone,1));

for ivar = 1:size(modelvar,1)
    imodelvar = modelvar(ivar,:);
    ivarname = varname{ivar};
    
    R_bftone = [];P_bftone = [];
    R_aftone = [];P_aftone = [];
    R_bfanswer = [];P_bfanswer = [];
    R_afanswer = [];P_afanswer = [];
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
    
%     figure;
%     scatter(R_bfanswer, R_afanswer);
%     xlabel('Coefficient with neural activity before answer');
%     ylabel('Coefficient with neural activity after answer');
%     title(['Pearson coefficient with ',ivarname]);
%     saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_answer.png']);
%     saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_answer.fig']);
%     close;
% 
%     figure;
%     scatter(R_bftone, R_aftone);
%     xlabel('Coefficient with neural activity before tone');
%     ylabel('Coefficient with neural activity after tone');
%     title(['Pearson coefficient with ',ivarname]);
%     saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_tone.png']);
%     saveas(gcf, [np_data_pathway, '\corr\scatter_',ivarname,'_tone.fig']);
%     close;

%     nunit = size(np2answer,1);
%     per01(ivar,1:4) = [sum(P_bftone<0.01) sum(P_aftone<0.01) sum(P_bfanswer<0.01) sum(P_afanswer<0.01)]/nunit;
%     per05(ivar,1:4) = [sum(P_bftone<0.05) sum(P_aftone<0.05) sum(P_bfanswer<0.05) sum(P_afanswer<0.05)]/nunit;

    
    %plot small time window
    ave_nbin =10;
    R2answer = [];
    R2tone = [];
    for iunit = 1:size(np2answer,1) 
        for jbin = 1:size(np2tone,3)/ave_nbin
            temp = squeeze(mean(np2tone(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
            [R2tone(iunit,jbin), P2tone(iunit,jbin)]= corr(temp', imodelvar','Type','Pearson');        
        end
        for jbin = 1:size(np2answer,3)/ave_nbin
            temp = squeeze(mean(np2answer(iunit,:,jbin*ave_nbin-9:jbin*ave_nbin), 3));
            [R2answer(iunit,jbin), P2answer(iunit,jbin)]= corr(temp', imodelvar','Type','Pearson');        
        end
    end

    R2tone_sorted = [];
    R2answer_sorted = [];
    R2tone(isnan(R2tone(:))) = 0;
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
    set(gca, 'Xtick', [0:10:size(np2answer,3)/ave_nbin], 'XTicklabel', ...
        100*[0:10:size(np2answer,3)/ave_nbin]);
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
    set(gca, 'Xtick', [0:10:2000/(bin_length*ave_nbin)], 'XTicklabel', ...
        [0:10*(bin_length*ave_nbin):2000]);
    title(['Pearson coefficient with ', ivarname,' 2tone'],'FontSize', 14);
    %tone time
    ph = patch([bin_stim_onset bin_stim_onset bin_stim_onset + 30 bin_stim_onset + 30]/ave_nbin, ...
          [size(R2tone, 1) 0 0  size(R2tone, 1)], 'w', 'edgecolor', 'none');
    alpha(ph,.4);
    filename = fullfile(np_data_pathway, 'corr_flAb', ivarname);
    saveas(gcf,[filename  '.png']);
    saveas(gcf,[filename  '.fig']);

end
   
R_allvar_bftone(isnan(R_allvar_bftone)) = 0;
R_allvar_aftone(isnan(R_allvar_aftone)) = 0;
R_allvar_bfanswer(isnan(R_allvar_bfanswer)) = 0;
R_allvar_afanswer(isnan(R_allvar_afanswer)) = 0;

save(fullfile(np_data_pathway, 'corr_flAb','corr.mat'),'varname',...
    'R_allvar_bftone','R_allvar_aftone','R_allvar_bfanswer','R_allvar_afanswer',...
     'P_allvar_bftone','P_allvar_aftone','P_allvar_bfanswer','P_allvar_afanswer');

%% plot scatter and fraction neuron with brain area
%ClusAbs_filtered(isnan(ClusAbs_filtered)) = max(ClusAbs_filtered(~isnan(ClusAbs_filtered)));
Abs = unique(ClusAbs_filtered);
Ab_idx = cell(1,length(Abs));
for iAb = 1:length(Abs)
    Ab_idx{iAb} = find(ClusAbs_filtered == iAb);
end

per01_areas = nan(length(Abs),4);%fraction of neurons significantly correlated in different areas
% plot
for jvar = 1:length(varname)
    figure;
    hold on
    for iAb = 1:length(Abs)
        per01_areas(iAb,3) = sum(P_allvar_bfanswer(jvar,Ab_idx{iAb})<0.01)/length(Ab_idx{iAb});
        per01_areas(iAb,4) = sum(P_allvar_afanswer(jvar,Ab_idx{iAb})<0.01)/length(Ab_idx{iAb});
        scatter(R_allvar_bfanswer(jvar,Ab_idx{iAb}),R_allvar_afanswer(jvar,Ab_idx{iAb}),20,'filled');    
    end
    limset = max(max(xlim, ylim));
    xlim([-1*limset limset]);
    ylim([-1*limset limset]);
    xlabel('Coefficient with neural activity before answer');
    ylabel('Coefficient with neural activity after answer');
    title(['Pearson coefficient with ',varname{jvar}]);

    xdiag = [-1*limset:0.05:limset];
    ydiag = xdiag;
    plot(xdiag,ydiag,'k--');

    legend(usedAb,'Location','Best');
    legend('boxoff');

    [r, p] = corr(R_allvar_bfanswer(jvar,:)',R_allvar_afanswer(jvar,:)');
    txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
    text(360,320, txtr, 'Units', 'pixels');
    text(360,300, txtp, 'Units', 'pixels');
    hold off;
    saveas(gcf, [np_data_pathway, '\corr_flAb\scatter_',varname{jvar},'_answer_area.png']);
    saveas(gcf, [np_data_pathway, '\corr_flAb\scatter_',varname{jvar},'_answer_area.fig']);
    close;

    %aline2tone
    figure;
    hold on
    for iAb = 1:length(Abs)
        per01_areas(iAb,1) = sum(P_allvar_bftone(jvar,Ab_idx{iAb})<0.01)/length(Ab_idx{iAb});
        per01_areas(iAb,2) = sum(P_allvar_aftone(jvar,Ab_idx{iAb})<0.01)/length(Ab_idx{iAb});
        scatter(R_allvar_bftone(jvar,Ab_idx{iAb}),R_allvar_aftone(jvar,Ab_idx{iAb}),20,'filled');    
    end
    limset = max(max(xlim, ylim));
    xlim([-1*limset limset]);
    ylim([-1*limset limset]);
    xlabel('Coefficient with neural activity before tone');
    ylabel('Coefficient with neural activity after tone');

    xdiag = [-1*limset:0.05:limset];
    ydiag = xdiag;
    plot(xdiag,ydiag,'k--');

    legend(usedAb,'Location','Best');
    legend('boxoff');

    [r, p] = corr(R_allvar_bftone(jvar,:)',R_allvar_aftone(jvar,:)');
    txtp = ['p = ', num2str(p,2)];txtr = ['r = ', num2str(r,2)];
    text(360,320, txtr, 'Units', 'pixels');
    text(360,300, txtp, 'Units', 'pixels');
    hold off;
    title(['Pearson coefficient with ',varname{jvar}]);

    saveas(gcf, [np_data_pathway, '\corr_flAb\scatter_',varname{jvar},'_tone_area.png']);
    saveas(gcf, [np_data_pathway, '\corr_flAb\scatter_',varname{jvar},'_tone_area.fig']);
    close;

    %show fraction of neurons significantly correlated
    figure;
    %set(gcf, 'Position', [100 100 700 600]);
    hold on
    x = 1:4;
    y = per01_areas;
    for i =1:length(Abs)
        plot(x,y(i,:),'.-','LineWidth',2, 'markersize',20);
    end
    legend(usedAb,'Location','best','FontSize',12);
    legend('boxoff');

    ylabel('Fraction of units','FontSize',20);
    title(['Fraction of units significantly correlated with ',varname{jvar}]);
    set(gca,'xTick',[1:4],'xTickLabel',{'bftone','aftone','bfanswer','afanswer'}, 'FontSize',10);

    saveas(gcf,fullfile(np_data_pathway, 'corr_flAb',  [varname{jvar} '_per01.png']));
    saveas(gcf,fullfile(np_data_pathway, 'corr_flAb',  [varname{jvar} '_per01.fig']));
end
%% plot summary for different brain areas
num_p01 = nan(length(varname),length(Abs),4);
for jvar = 1:length(varname)
    for iAb = 1:length(Abs)
         num_p01(jvar, iAb, 1) = sum(P_allvar_bftone(jvar,Ab_idx{iAb})<0.01);
         num_p01(jvar, iAb, 2) = sum(P_allvar_aftone(jvar,Ab_idx{iAb})<0.01);
         num_p01(jvar, iAb, 3) = sum(P_allvar_bfanswer(jvar,Ab_idx{iAb})<0.01);
         num_p01(jvar, iAb, 4) = sum(P_allvar_afanswer(jvar,Ab_idx{iAb})<0.01);
    end
end

var_interest = [1;2;5;6];%delta, P_boundary low(t+1), P_choiceleft, P_boundary high(t+1)
varname_interest = {varname{1}, varname{2}, varname{5}, varname{6}};
for iAb = 1:length(Abs)
    
    figure;
    hold on
    x = 1:4;
    y = squeeze(num_p01(var_interest, iAb,:));
    for i =1:length(var_interest)
        plot(x,y(i,:),'.-','LineWidth',2, 'markersize',20);
    end
    legend(varname_interest,'Location','best','FontSize',12);
    legend('boxoff');

    ylabel('Number of units','FontSize',20);
    title(usedAb{iAb});
    set(gca,'xTick',[1:4],'xTickLabel',{'bftone','aftone','bfanswer','afanswer'}, 'FontSize',10);

    saveas(gcf,fullfile(np_data_pathway, 'corr_flAb',  [usedAb{iAb} '_per01.png']));
    saveas(gcf,fullfile(np_data_pathway, 'corr_flAb',  [usedAb{iAb} '_per01.fig']));
end


%% plot example neuron
iinterest = 2;
x6 = find(diff(BlockType(inds_use)));
varinterest = modelvar(iinterest,:);
[~,example_idx] = sort(P_allvar_afanswer(iinterest,:));
for ie = 1:5
    
    figure;
    ie_idx = example_idx(ie);    
    example_ID = ClusIDs_filtered(ie_idx);
    
    sgtitle(['unit ' num2str(example_ID) ' ' usedAb{ClusAbs_filtered(ie_idx)}]);
    set(gcf, 'Position', [100 100 1600 700]);
    subplot(1,2,1)
    x = 1:length(varinterest);
    y1 = smooth(zscore(varinterest),10);
    %     y2 = smooth(zscore(mean(squeeze(np2tone(example_num,:,1:bin_stim_onset-1)),2)),10);
    %     y3 = smooth(zscore(mean(squeeze(np2tone(example_num,:,bin_stim_onset:200)),2)),10);
    %     y4 = smooth(zscore(mean(squeeze(np2answer(example_num,:,1:bin_answer-1)),2)),10);
    y5 = smooth(zscore(mean(squeeze(np2answer(ie_idx,:,bin_answer:end)),2)),10);
    plot(x,y1,'k','LineWidth',1);
    hold on;
    %     plot(x,y2,'y','LineWidth',1.5);
    %     plot(x,y3,'g','LineWidth',1.5);
    %     plot(x,y4,'c','LineWidth',1.5);
    plot(x,y5,'b','LineWidth',1);
    
    yl = ylim;
    y6 = linspace(yl(1),yl(2), 20);
    for iblock_change = 1:length(x6)
        plot(x6(iblock_change)*ones(20,1),y6,'k--','LineWidth',1);% block change
    end
    
    xlabel('# Trial','FontSize',12);
    ylabel('Z score','FontSize',12);
    legend('Model','Neural activity afanswer','FontSize',12, 'Location', 'best');
    legend('boxoff');
    hold off;

    subplot(1,2,2);
    scatter(y1,y5);
    xlabel(['z-scored ',varname(iinterest)],'FontSize',12);
    ylabel('z-scored firing rate','FontSize',12);
    saveas(gcf,fullfile(np_data_pathway, 'corr_flAb', [varname{iinterest}, ' unit ', num2str(example_ID),'.png']));
    saveas(gcf,fullfile(np_data_pathway, 'corr_flAb', [varname{iinterest}, ' unit ', num2str(example_ID),'.fig']));
    close;
end 
 
%  %% plot venn
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