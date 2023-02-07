function UsedTrigOnTime = triggerOnsetTime(ksFolder,TrigwaveScales,trigTime_ms,sample_rate)
% TriggerType should be a string variable indicates what the
% trigger onset time is indicated, either task or passive
% sessions
% used method: UsedTrigOnTime = triggerOnsetTime(ksFolder,[],[6,2],30000);
% % trigTime_ms should be [4] for passive sessions

% trigTimeLen should in ms form
if isempty(trigTime_ms)
    trigTime_ms = 2; % default is 2ms
end

if isempty(TrigwaveScales)
    try
        trigScaleStrc = load(fullfile(ksFolder,'..','TriggerDatas.mat'),'TriggerEvents');
    catch
        trigScaleStrc = load(fullfile(ksFolder,'TriggerDatas.mat'),'TriggerEvents');
    end
    TrigwaveScales = trigScaleStrc.TriggerEvents;
end

trigLen = trigTime_ms/1000*sample_rate;
TrigWaveAll_lens = TrigwaveScales(:,2) - TrigwaveScales(:,1);
TrigWaveAll_lens(TrigWaveAll_lens < 0.0005*sample_rate) = []; % exclude some unknowed zeros trigger events

if length(trigLen) > 1
    % for task condition, the first 10 trial will have different trigger durations
    Trig_InitSeg_durs = (abs(TrigWaveAll_lens - trigLen(1)) < 5);
    Trig_InitSeg_consec = consecTRUECount(Trig_InitSeg_durs);
    Trig_Init_TrialEnd = find(Trig_InitSeg_consec == 10,1,'first'); % used the first 10 tirals
    if isempty(Trig_Init_TrialEnd)
        warning('Failed to find enough lengthed consecutive trials.\n Please check your trigger duration data or input value.\n')
    end
    SessionTrStartInds = find(abs(TrigWaveAll_lens((Trig_Init_TrialEnd(1)+1):end) - trigLen(2)) < 5 | ...
        abs(TrigWaveAll_lens((Trig_Init_TrialEnd(1)+1):end) - trigLen(1)) < 5)+Trig_Init_TrialEnd(1);
    TrigWaveAll_lensEquals = [((Trig_Init_TrialEnd-9):Trig_Init_TrialEnd)';SessionTrStartInds];
    
    UsedTrigOnTime = TrigwaveScales(TrigWaveAll_lensEquals,1)/sample_rate; % in seconds
else
    % for passive condition, unique trigger duration
    TrigWaveAll_lensEquals = abs(TrigWaveAll_lens - trigLen) < 2;
    UsedTrigOnTime= TrigwaveScales(TrigWaveAll_lensEquals,1)/sample_rate; % in seconds
end


fprintf('Totally %d number of triggers were detected.\n',length(UsedTrigOnTime));

end
