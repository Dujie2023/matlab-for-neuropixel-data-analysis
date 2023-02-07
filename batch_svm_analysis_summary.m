%after svm analysis for each probe, organiza the results according to
%regions

results_savefolder = 'D:\DJ_data\catprocessdata_bin_path\auc_analysis_summary';
%batch data pathway
filepath_path = 'D:\DJ_data\catprocessdata_bin_path';
filepath_name = 'npdata_process_path_202206batch.xlsx';
[~,filepathtemp,~] = xlsread(fullfile(filepath_path,filepath_name));
filenum = length(filepathtemp)-1;

%output data
load('AbstractBrainRegions.mat');
SVMSummaryRegions = abinfo;
% field: 1, acc for all units(nsessions*4 matrix, 4 time windows); 
%2. acc for 10 units (nsessions*4 matrix, 4 time windows);
%3. acc for small timewindows,alined 2 tone(nsessions*ntime windows);
%4. acc for small timewindows,alined 2 answer(nsessions*ntime windows);
%5-8.similar to 1-4, but for different labels
%9. unit number(nsessios*1);

abareanum = size(SVMSummaryRegions,1);
SVMSummaryRegions(1).accblock = [];
SVMSummaryRegions(1).acc10block = [];
SVMSummaryRegions(1).acc100ms2toneblock = [];
SVMSummaryRegions(1).acc100ms2answerblock = [];
SVMSummaryRegions(1).accchoice = [];
SVMSummaryRegions(1).acc10choice = [];
SVMSummaryRegions(1).acc100ms2tonechoice = [];
SVMSummaryRegions(1).acc100ms2answerchoice = [];
SVMSummaryRegions(1).unitnum = [];
%% organize data from all files
for ifile = 1:filenum
    
    np_data_pathway = fullfile(filepathtemp{ifile+1},'ks2_5');
    %np_data_pathway =  'D:\DJ_data\catprocessdata\202206batch_mouse1_cat\Cat_20220828_mouse1_OFC1R_OFC2L_PPC2L_g0_imec0\ks2_5';
    areainfo_filename = 'probeAbstractBrainRegions.mat';
    load(fullfile(np_data_pathway,areainfo_filename));
    
    try
       aucinfo_filename = 'BlockSVM_acc';  
        load(fullfile(np_data_pathway,'svm',aucinfo_filename));    
        BlockAcc = probeAreaAcc;
        aucinfo_filename = 'ChoiceSVM_acc';
        load(fullfile(np_data_pathway,'svm',aucinfo_filename));    
        ChoiceAcc = probeAreaAcc;
    catch
        fprintf("svm analysis results don't exist: \n");
        disp(np_data_pathway);
        continue;
    end
    
    if ~isfield(probeAreaAcc,'accbftone')
        fprintf('Unit number is not enough: \n');
        disp(np_data_pathway);
        continue        
    end
                
    for i = 1:abareanum

        accblock = [BlockAcc(i).accbftone BlockAcc(i).accaftone...
            BlockAcc(i).accbfanswer BlockAcc(i).accafanswer];   
        acc10block = [BlockAcc(i).acc10bftone BlockAcc(i).acc10aftone...
            BlockAcc(i).acc10bfanswer BlockAcc(i).acc10afanswer];
        acc100ms2toneblock = BlockAcc(i).acc100ms2tone;
        acc100ms2answerblock = BlockAcc(i).acc100ms2answer;    

        accchoice = [ChoiceAcc(i).accbftone ChoiceAcc(i).accaftone...
            ChoiceAcc(i).accbfanswer ChoiceAcc(i).accafanswer];
        acc10choice = [ChoiceAcc(i).acc10bftone ChoiceAcc(i).acc10aftone...
            ChoiceAcc(i).acc10bfanswer ChoiceAcc(i).acc10afanswer];
        acc100ms2tonechoice = ChoiceAcc(i).acc100ms2tone;
        acc100ms2answerchoice = ChoiceAcc(i).acc100ms2answer;

    
        SVMSummaryRegions(i).accblock = [SVMSummaryRegions(i).accblock; ...
            accblock];
        SVMSummaryRegions(i).acc10block = [SVMSummaryRegions(i).acc10block;...
            acc10block];
        SVMSummaryRegions(i).acc100ms2toneblock = [SVMSummaryRegions(i).acc100ms2toneblock;...
            acc100ms2toneblock];
        SVMSummaryRegions(i).acc100ms2answerblock = [SVMSummaryRegions(i).acc100ms2answerblock;...
            acc100ms2answerblock];
        SVMSummaryRegions(i).accchoice = [SVMSummaryRegions(i).accchoice;...
            accchoice];
        SVMSummaryRegions(i).acc10choice = [SVMSummaryRegions(i).acc10choice;...
            acc10choice];
        SVMSummaryRegions(i).acc100ms2tonechoice = [SVMSummaryRegions(i).acc100ms2tonechoice;...
            acc100ms2tonechoice];
        SVMSummaryRegions(i).acc100ms2answerchoice = [SVMSummaryRegions(i).acc100ms2answerchoice;...
            acc100ms2answerchoice];
        
        if probeAreaAcc(i).probeUsedClus_num>=10
            SVMSummaryRegions(i).unitnum = [SVMSummaryRegions(i).unitnum;...
                probeAreaAcc(i).probeUsedClus_num];
        end
    end 
end

%% plot
%figure 1: plot for each area
for i  = 1:abareanum
    if ~isempty(SVMSummaryRegions(i).accblock)
        h = figure;
        titletext = [SVMSummaryRegions(i).safe_name,'  (nsession=',num2str(length(SVMSummaryRegions(i).unitnum)),')'];
        sgtitle(titletext);

        try
        %boxplot: acc for Block and Choice
        subplot(2,1,1);        
        boxplot(SVMSummaryRegions(i).accblock,'Labels',{'bf tone','af tone','bf answer','af answer'});
        ylim([0.5 1]);
        ylabel('Decoding Accuracy for Block','FontSize',10);
        box off;
        
        subplot(2,1,2);
        boxplot(SVMSummaryRegions(i).accchoice,'Labels',{'bf tone','af tone','bf answer','af answer'});
        ylim([0.5 1]);
        ylabel('Decoding Accuracy for Choice','FontSize',10);
        box off;
        
        end
        saveas(gcf,fullfile(results_savefolder,'regions_svm',[num2str(i),'.png']));        
        saveas(gcf,fullfile(results_savefolder,'regions_svm',[num2str(i),'.fig']));        
        saveppt(fullfile(results_savefolder,'svm.ppt'),' ',h)
        close;
        
        
        %plot: acc in small time windows for Block
        %plot: acc in small time windows for Choice
    end
end

%figure2: summary 
%auc for all recorded regions
