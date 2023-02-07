%after auc analysis for each probe, organiza the results according to
%regions

results_savefolder = 'D:\DJ_data\catprocessdata_bin_path\auc_analysis_summary';
%batch data pathway
filepath_path = 'D:\DJ_data\catprocessdata_bin_path';
filepath_name = 'npdata_process_path_202206batch.xlsx';
[~,filepathtemp,~] = xlsread(fullfile(filepath_path,filepath_name));
filenum = length(filepathtemp)-1;

%output data
load('AbstractBrainRegions.mat');
AUCSummaryRegions = abinfo;% field: 1, auc(nunits*4 matrix, 4 time windows); 
%2. corresponding sig 99% threshold(nunits*4 matrix);3-6: auc and sig
%threshold, but for diffferent labels;
%7. use trial number (nunits*1); 8. switch times(nunit*1);
abareanum = size(AUCSummaryRegions,1);
AUCSummaryRegions(1).aucblock = [];
AUCSummaryRegions(1).aucblock_sigthre = [];
AUCSummaryRegions(1).aucchoice = [];
AUCSummaryRegions(1).aucchoice_sigthre = [];
AUCSummaryRegions(1).aucoutcome = [];
AUCSummaryRegions(1).aucoutcome_sigthre = [];
AUCSummaryRegions(1).trinum = [];
AUCSummaryRegions(1).switchtimes = [];

%% session data loading
for ifile = 1:filenum
    
    np_data_pathway = fullfile(filepathtemp{ifile+1},'ks2_5');
    %np_data_pathway =  'D:\DJ_data\catprocessdata\202206batch_mouse1_cat\Cat_20220828_mouse1_OFC1R_OFC2L_PPC2L_g0_imec0\ks2_5';
    areainfo_filename = 'probeAbstractBrainRegions.mat';
    load(fullfile(np_data_pathway,areainfo_filename));

    %block infomation
    aucinfo_filename = 'Block auc';
    try
        load(fullfile(np_data_pathway,'auc',aucinfo_filename));    
    catch
        fprintf("auc analysis results don't exist: \n");
        disp(np_data_pathway);
        continue;
    end
    
    %auc data organization
    auc_block = [auc_bftone auc_aftone auc_bfanswer auc_afanswer];
    auc_blocksigthre = [sig_thre01_aline2tone sig_thre01_aline2tone sig_thre01_aline2answer sig_thre01_aline2answer];
    trinumber = length(cur_label);
    switchtimes = sum(diff(cur_label)~=0);

    %choice infomation
    aucinfo_filename = 'Choice auc';
    load(fullfile(np_data_pathway,'auc',aucinfo_filename));
    choice_auc = [auc_bftone auc_aftone auc_bfanswer auc_afanswer];
    auc_choicesigthre = [sig_thre01_aline2tone sig_thre01_aline2tone sig_thre01_aline2answer sig_thre01_aline2answer];

    %outcome infomation
    aucinfo_filename = 'Outcome auc';
    load(fullfile(np_data_pathway,'auc',aucinfo_filename));
    outcome_auc = [auc_bftone auc_aftone auc_bfanswer auc_afanswer];
    auc_outcomesigthre = [sig_thre01_aline2tone sig_thre01_aline2tone sig_thre01_aline2answer sig_thre01_aline2answer];

    for i = 1:abareanum
        AUCSummaryRegions(i).aucblock = [AUCSummaryRegions(i).aucblock;...
            auc_block(probeAbArea(i).idx_probeUsedClus_IDs,:)];
        AUCSummaryRegions(i).aucblock_sigthre = [AUCSummaryRegions(i).aucblock_sigthre; ...
            ones(probeAbArea(i).probeUsedClus_num,1)*auc_blocksigthre];

        AUCSummaryRegions(i).aucchoice = [AUCSummaryRegions(i).aucchoice;...
            choice_auc(probeAbArea(i).idx_probeUsedClus_IDs,:)];
        AUCSummaryRegions(i).aucchoice_sigthre = [AUCSummaryRegions(i).aucchoice_sigthre; ...
            ones(probeAbArea(i).probeUsedClus_num,1)*auc_choicesigthre];

        AUCSummaryRegions(i).aucoutcome = [AUCSummaryRegions(i).aucoutcome;...
            outcome_auc(probeAbArea(i).idx_probeUsedClus_IDs,:)];
        AUCSummaryRegions(i).aucoutcome_sigthre = [AUCSummaryRegions(i).aucoutcome_sigthre; ...
            ones(probeAbArea(i).probeUsedClus_num,1)*auc_outcomesigthre];

        AUCSummaryRegions(i).trinum = [AUCSummaryRegions(i).trinum; ...
            ones(probeAbArea(i).probeUsedClus_num,1)*trinumber];
        AUCSummaryRegions(i).switchtimes = [AUCSummaryRegions(i).switchtimes; ...
            ones(probeAbArea(i).probeUsedClus_num,1)*switchtimes];
    end

end

%% plot for different regions
%figure 1
for i = 1:abareanum
    if length(AUCSummaryRegions(i).trinum)>=20 %unit number
        h = figure;
        set(gcf,'Position',[100,50,900,900]);
        titletext = [AUCSummaryRegions(i).safe_name,'  (n=',num2str(length(AUCSummaryRegions(i).trinum)),')'];
        sgtitle(titletext);
   
        %bftone block vs aftone block auc
        subplot(3,3,1);
        aucblock = AUCSummaryRegions(i).aucblock;
        blockthre = AUCSummaryRegions(i).aucblock_sigthre;
        sigidx = {aucblock(:,1)>blockthre(:,1),aucblock(:,2)>blockthre(:,1)};
        gscatter(aucblock(:,1),aucblock(:,2),sigidx,'kgbr','oooo');
        xlabel('AUC bf tone for Block');
        ylabel('AUC af tone for Block');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');
        box off;
        legend('hide');
        
        %bfanswer block auc vs afanswer block auc
        subplot(3,3,2);
        sigidx = {aucblock(:,3)>blockthre(:,3),aucblock(:,4)>blockthre(:,3)};
        gscatter(aucblock(:,3),aucblock(:,4),sigidx,'kgbr','oooo');
        xlabel('AUC bf answer for Block');
        ylabel('AUC af answer for Block');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');
        box off;
        legend('hide');
        
        %bftone choice auc vs aftone choice auc
        subplot(3,3,4);
        aucchoice = AUCSummaryRegions(i).aucchoice;
        choicethre = AUCSummaryRegions(i).aucchoice_sigthre;
        sigidx = {aucchoice(:,1)>choicethre(:,1),aucchoice(:,2)>choicethre(:,1)};
        gscatter(aucchoice(:,1),aucchoice(:,2),sigidx,'kgbr','oooo');
        xlabel('AUC bf tone for Choice');
        ylabel('AUC af tone for Choice');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');
        box off;
        legend('hide');
        
        %bfanswer choice AUC vs afanswer choice AUC
        subplot(3,3,5);
        sigidx = {aucchoice(:,3)>choicethre(:,3),aucchoice(:,4)>choicethre(:,3)};
        gscatter(aucchoice(:,3),aucchoice(:,4),sigidx,'kgbr','oooo');
        xlabel('AUC bf answer for Choice');
        ylabel('AUC af answer for Choice');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');
        box off;
        legend('hide');
        
        %bftone outcome AUC vs aftone outcome AUC
        subplot(3,3,7);
        aucoutcome = AUCSummaryRegions(i).aucoutcome;
        outcomethre = AUCSummaryRegions(i).aucoutcome_sigthre;
        sigidx = {aucoutcome(:,1)>outcomethre(:,1),aucoutcome(:,2)>outcomethre(:,2)};
        gscatter(aucoutcome(:,1),aucoutcome(:,2),sigidx,'kgbr','oooo');
        xlabel('AUC bf tone for Outcome');
        ylabel('AUC af tone for Outcome');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');
        box off;
        legend('hide');
        
        %bfanswer outcome AUC vs afanswer outcome AUC
        subplot(3,3,8);
        sigidx = {aucoutcome(:,3)>outcomethre(:,3),aucoutcome(:,4)>outcomethre(:,3)};
        gscatter(aucoutcome(:,3),aucoutcome(:,4),sigidx,'kgbr','oooo');
        xlabel('AUC bf answer for Outcome');
        ylabel('AUC af answer for Outcome');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');
        box off;
        legend('hide');
        
        %bftone block AUC vs aftone choice
        subplot(3,3,3);
        sigidx = {aucblock(:,1)>blockthre(:,1),aucchoice(:,2)>choicethre(:,1)};        
        gscatter(aucblock(:,1),aucchoice(:,2),sigidx,'kgbr','oooo');
        xlabel('AUC bf tone for Block');
        ylabel('AUC af tone for Choice');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');    
        box off;
        legend('hide');
        
        %afanswer outcome AUC vs afanswer choice AUC
        subplot(3,3,6);
        sigidx = {aucchoice(:,4)>choicethre(:,3),aucoutcome(:,4)>outcomethre(:,3)};        
        gscatter(aucchoice(:,4),aucoutcome(:,4),sigidx,'kgbr','oooo');
        xlabel('AUC af answer for Choice');
        ylabel('AUC af answer for Outcome');
        xr = xlim;yr = ylim;
        newrange= [min(xr(1),yr(1)),max(xr(2),yr(2))];
        xlim(newrange);ylim(newrange);
        set(gca,'TickDir','out');    
        box off;
        legend('hide');
        
        %switch time histogram
        subplot(3,3,9);
        histogram(AUCSummaryRegions(i).switchtimes,'FaceColor','k','EdgeColor','w');
        xlabel('Switchtimes');
        ylabel('#units');
        set(gca,'XLim',[0.5 5.5]);
        set(gca,'XTick',[1,2,3,4,5]);
        set(gca,'TickDir','out');
        box off
        saveas(gcf,fullfile(results_savefolder,'regions_auc',[num2str(i),'.png']));        
        saveas(gcf,fullfile(results_savefolder,'regions_auc',[num2str(i),'.fig']));        
        saveppt(fullfile(results_savefolder,'auc.ppt'),' ',h)
        close;

    end
end


%% 
% unitnum = 0;
% for i = 1:105
% unitnum=unitnum+length(AUCSummaryRegions(i).trinum);
% end