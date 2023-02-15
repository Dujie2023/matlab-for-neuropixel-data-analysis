%auc analysis for neuropixel data
%cal cuc for 3 labels: Chocie, Block, Outcome
%in 4 time windows:1s befor tone stim onset, 1s after tone stim onset, 1s
%before action, 2s after action
filepath_path = 'D:\DJ_data\catprocessdata_bin_path';
filepath_name = 'npdata_process_path_202206batch.xlsx';
[~,filepathtemp,~] = xlsread(fullfile(filepath_path,filepath_name));
filenum = length(filepathtemp)-1;

for ifile = 1:filenum
    %% loading data
    %np_data_pathway =  'D:\DJ_data\catprocessdata\202206batch_mouse1_cat\20220828_mouse1_OFC1R_OFC2L_PPC2L_g0_cat\Cat_20220828_mouse1_OFC1R_OFC2L_PPC2L_g0_imec0\ks2_5';
    np_data_pathway = fullfile(filepathtemp{ifile+1},'ks2_5');
    
    fprintf('processing %d/%d file...\n',ifile,filenum);
    preAUCprocessing_DJ;
    
    load(fullfile(np_data_pathway,'BehavData.mat'));
    inds_use = behavResults.Action_choice~=2;
    blocktype_use = behavResults.BlockType(inds_use);
    choice_use = behavResults.Action_choice(inds_use);
    outcome = behavResults.Time_reward~=0;
    outcome_use = outcome(inds_use);

    %% auc analysis
    nbin_aftone = 100;%1s after tone stim onset
    nbin_afanswer = 200;%2s after action;

%     save_folder = fullfile(np_data_pathway, 'auc');    
%     if exist(save_folder)
%         recycle('off');        
%         dinfo = dir(save_folder);
%         dinfo([dinfo.isdir]) = [];   %skip directories
%         filenames = fullfile(save_folder, {dinfo.name});
%         if ~isempty(filenames)
%             delete(filenames{:})
%         end
%     end    
    
    auc_session(np_data_pathway,blocktype_use,'Block',nbin_aftone, nbin_afanswer,1);   
    auc_session(np_data_pathway,choice_use,'Choice',nbin_aftone, nbin_afanswer,1);   
    auc_session(np_data_pathway,outcome_use,'Outcome',nbin_aftone, nbin_afanswer,1);

end


