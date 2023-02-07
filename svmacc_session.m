function svmacc_session(np_data_pathway, cur_label,labelname,nbin_aftone, nbin_afanswer)   
%analyze neural acitivity recorded by neuropixels using SVM. 
%the function processes recorded neural activity data in a session, and
%cal the accuracy of SVM classifier for each brain area. Brain regions in
%which unit number<10 will not be used.
%analysis results will be saved in a struct(nareas*  fields)

%-----------------input----------------
%np_data_pathway: ks2_5 folder pathway
%cur_label: label vector(ntrials), e.g Choice, Outcome, Block
%label name: char
%------time window setting
%nbin_aftone: how many bins after stime tone onset are used, 10ms for a bin;
%nbin_afanswer: how many bins after answer are used, 10ms for a bin;

%% data loading
load(fullfile(np_data_pathway, 'np_extracted.mat'));

%if unit number is 0, skip analysis
if isempty(Spike_firingrate_bin)
    warning('The good unit number of this probe is 0.');
end

save_folder = fullfile(np_data_pathway, 'svm');
if ~exist(save_folder)
    mkdir(save_folder);
end
savefilename = [labelname 'SVM_acc.mat'];

np2tone = Spike_firingrate_bin;%neural activity, alined to tone onset time, matrix(nunits*ntrials*nbins),unit:Hz
np2answer = Spike_firingrate_bin_outcome; %neural activity, alined to answer time, matrix(nunits*ntrials*nbins)
neuroactivity_bftone = mean(np2tone(:,:,1:bin_stim_onset-1), 3);
neuroactivity_aftone = mean(np2tone(:,:,bin_stim_onset:bin_stim_onset+nbin_aftone), 3);    
neuroactivity_bfanswer = mean(np2answer(:,:,1:bin_answer-1), 3);
neuroactivity_afanswer = mean(np2answer(:,:,bin_answer:bin_answer+nbin_afanswer), 3);    

%% load area infomation
areainfo_filename = 'probeAbstractBrainRegions.mat';
load(fullfile(np_data_pathway,areainfo_filename));

minUnitNumber = 10;%Brain regions in which unit number<10 will not be used.
%% cal SVM accuracy

probeAreaAcc = probeAbArea;

for iArea = 1:size(probeAbArea,1)
    area_unit_idx = probeAbArea(iArea).idx_probeUsedClus_IDs;

    if length(area_unit_idx)>=minUnitNumber
        probeAreaAcc(iArea).accbftone = cal_svmacc_npdata(neuroactivity_bftone(area_unit_idx,:)', cur_label');
        probeAreaAcc(iArea).accaftone = cal_svmacc_npdata(neuroactivity_aftone(area_unit_idx,:)', cur_label');
        probeAreaAcc(iArea).accbfanswer = cal_svmacc_npdata(neuroactivity_bfanswer(area_unit_idx,:)', cur_label');
        probeAreaAcc(iArea).accafanswer = cal_svmacc_npdata(neuroactivity_afanswer(area_unit_idx,:)', cur_label');

        probeAreaAcc(iArea).acc10bftone = calmeanSVMAcc(neuroactivity_bftone(area_unit_idx,:)', cur_label');
        probeAreaAcc(iArea).acc10aftone = calmeanSVMAcc(neuroactivity_aftone(area_unit_idx,:)', cur_label');
        probeAreaAcc(iArea).acc10bfanswer = calmeanSVMAcc(neuroactivity_bfanswer(area_unit_idx,:)', cur_label');
        probeAreaAcc(iArea).acc10afanswer = calmeanSVMAcc(neuroactivity_afanswer(area_unit_idx,:)', cur_label');

    end
end

%% cal acc for small timewindow(10bins, 100ms)

t_bftone = 1;%s
t_aftone = 3;%s

t_bfanswer = 1;
t_afanswer = 2;
numbin  = 10;% small time window length: numbin*bin_length = 10*10ms = 100ms
idx_twpoint = 1:numbin:1000/bin_length*(t_bftone+t_aftone);
idx_twpoint2 = 1:numbin:1000/bin_length*(t_bfanswer+t_afanswer);

np2tonebin100 = nan(size(np2tone,1),size(np2tone,2),length(idx_twpoint));
for j = 1:length(idx_twpoint)
    np2tonebin100(:,:,j) = mean(np2tone(:,:,idx_twpoint(j):idx_twpoint(j)+9),3);
end

np2answerbin100 = nan(size(np2answer,1),size(np2answer,2),length(idx_twpoint2));
for j = 1:length(idx_twpoint2)
    np2answerbin100(:,:,j) = mean(np2answer(:,:,idx_twpoint2(j):idx_twpoint2(j)+9),3);
end

for iArea = 1:size(probeAbArea,1)
    area_unit_idx = probeAbArea(iArea).idx_probeUsedClus_IDs;
    if length(area_unit_idx)>=minUnitNumber
        acc = nan(1,length(idx_twpoint));
        acc2 = nan(1,length(idx_twpoint2));
        
        for t = 1:length(idx_twpoint)
           na = squeeze(np2tonebin100(area_unit_idx,:,t));
           acc(t)  = cal_svmacc_npdata(na', cur_label');
        end 
       
        for t = 1:length(idx_twpoint2)
           na = squeeze(np2answerbin100(area_unit_idx,:,t));
           acc2(t)  = cal_svmacc_npdata(na', cur_label');
        end 
        probeAreaAcc(iArea).acc100ms2tone = acc;    
        probeAreaAcc(iArea).acc100ms2answer = acc2;    

    end
    
end
%% save results
 save(fullfile(save_folder,savefilename), 'probeAreaAcc');
end

