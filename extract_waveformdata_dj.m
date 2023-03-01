function [] = extract_waveformdata_dj(np_data_pathway, bin_name, rez)
%function to extract waveform data from bin and kilosort output
%output:UnitDatas(num_clusters,2)
%date:2022/3/23

ksfolder = fullfile(np_data_pathway,'ks2_5');
if ~isfolder(ksfolder)
    ksfolder = fullfile(rez.ops.ksFolderPath,'kilosort3');
end

bin_path = fullfile(np_data_pathway,bin_name);
ftempid = fopen(bin_path);

if ~isfolder(fullfile(ksfolder,'Waveforms'))
    mkdir(fullfile(ksfolder,'Waveforms'));
end

% folders for saving process data
ksprocessfolder = fullfile(np_data_pathway,'ks2_5process');
if ~isfolder(ksprocessfolder)
    mkdir(fullfile(np_data_pathway,'ks2_5process');
end
if ~isfolder(fullfile(ksprocessfolder,'Waveforms'))
    mkdir(fullfile(ksprocessfolder,'Waveforms'));
end

% load cluster info
cgsFile = fullfile(ksfolder,'cluster_info.csv');

% cluster_info = readtable(cgsFile);
UsedClus_IDs = cluster_info.id;
ChannelUseds_id = cluster_info.ch;
SpikeTimeSample = readNPY(fullfile(ksfolder,'spike_times.npy'));
SpikeClus = readNPY(fullfile(ksfolder,'spike_clusters.npy'));

WaveWinSamples = [-30,51];%wafeform duration
NumofUnit = length(UsedClus_IDs);
UnitDatas = cell(NumofUnit,2);
UnitFeatures = cell(NumofUnit,7);

for cUnit = 1: NumofUnit
    
    cClusInds = UsedClus_IDs(cUnit);
    %cClusChannel = ChannelUseds_id(cUnit);
    cClusChannel = ChannelUseds_id(cUnit)+1;
    cClus_Sptimes = SpikeTimeSample(SpikeClus == cClusInds);
    
    %maximum of SPNums is 2000
    if numel(cClus_Sptimes) < 2000
        UsedSptimes = cClus_Sptimes;
        SPNums = length(UsedSptimes);
    else
        UsedSptimes = cClus_Sptimes(randsample(numel(cClus_Sptimes),2000));
        SPNums = 2000;
    end  
    
    cspWaveform = nan(SPNums,diff(WaveWinSamples));  %create matrix(SPNums, WaveWinSamples) with elements NaN
    AllChannelWaveData = nan(SPNums,384,diff(WaveWinSamples));
    for csp = 1 : SPNums
        cspTime = UsedSptimes(csp);
        cspStartInds = cspTime+WaveWinSamples(1);
        cspEndInds = cspTime+WaveWinSamples(2);
        offsetTimeSample = cspStartInds - 1;
        
        %ensure that the spike in the recording time
        if offsetTimeSample < 0 || cspEndInds > rez.ops.sampsToRead
            continue;
        end  
        
        %read waveform data from .bin file
        offsets = 385*(cspStartInds-1)*2;   %position indicator
        status = fseek(ftempid,offsets,'bof');
        if ~status
            % correct offset value is set
            AllChnDatas = fread(ftempid,[385 diff(WaveWinSamples)],'int16');% read data with dimensions [385 diff(WaveWinSamples)]
            cspWaveform(csp,:) = AllChnDatas(cClusChannel,:);% waveform data in a channel
            AllChannelWaveData(csp,:,:) = AllChnDatas(1:384,:); % for waveform spread calculation
        end
    end
    
    %save waveform data
    if size(cspWaveform,1) == 1
        AvgWaves = cspWaveform;
        UnitDatas{cUnit,2} = squeeze(AllChannelWaveData);
    else
        AvgWaves = mean(cspWaveform,'omitnan');%omitnan means the mean is of all non-NaN elements
        UnitDatas{cUnit,2} = squeeze(mean(AllChannelWaveData,'omitnan'));
    end
    UnitDatas{cUnit,1} = cspWaveform;
    
    clearvars AllChannelWaveData
    %% plot and calculate waveform feature
    huf = figure('visible','off');
    
    waveform_feature_funcollect = waveform_funcollect;
    [wave_amplitude, trough2peakT, repolarT, prepeak2pospeakratio, troughpeakidx] = ...
                                waveform_feature_funcollect.cal_waveform_feature(AvgWaves,WaveWinSamples);

    plot(AvgWaves);
    title(num2str(UsedClus_IDs(cUnit),'cluster=%d'));
    text(6,0.8*max(AvgWaves),{sprintf('trough2peak = %d',trough2peakT);...
        sprintf('repolarT = %d',repolarT);...
        sprintf('prepeak2pospeakratio = %.2f',prepeak2pospeakratio)},'FontSize',10);    
      
    saveName = fullfile(ksprocessfolder,'Waveforms_bin',sprintf('Unit%d waveform plot save',cUnit));
    saveas(huf,saveName,'png');   
    saveas(huf,saveName,'fig');   
    close(huf);
    
    fprintf('Unit%d/%d waveform plot save.\n',cUnit, NumofUnit)
end

save(fullfile(ksprocessfolder,'WaveformDatas.mat'), 'UnitDatas', 'UnitFeatures', '-v7.3');
fclose(ftempid);




