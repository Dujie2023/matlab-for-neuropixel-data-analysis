function [Time_alined_outcome, Time_alined_tone,aline_time_window, aline_time_window_outcome,...
    t_after_tone_onset, t_before_tone_onset, t_after_answer, t_before_answer] = ... 
    alinetime(Tone_onset_time, Answer_time,inds_use,TaskTriggerStartTime)

%aline time
t_before_tone_onset = 1000;%ms
t_after_tone_onset = 5000;%ms
Tone_onsetT_omiss = Tone_onset_time(inds_use);
MinOnsetTime = min(Tone_onsetT_omiss);
if MinOnsetTime <= t_before_tone_onset
    t_before_tone_onset = double(MinOnsetTime);
end

trial_all_time = 1000*diff(TaskTriggerStartTime);%unit:ms
trial_all_time(length(trial_all_time)+1,1) = 1000*10;

min_tone2endT = round(min(trial_all_time(inds_use)' - double(Tone_onsetT_omiss)));
if t_after_tone_onset >= min_tone2endT
    t_after_tone_onset = min_tone2endT;
end

for triN = 1:length(Tone_onsetT_omiss)
    Time_alined_tone(triN,:) = double(Tone_onsetT_omiss(triN)) + [-t_before_tone_onset+1 t_after_tone_onset];% time 1s before tone onset, time 6-8s after tone offset;
end
aline_time_window =  t_after_tone_onset + t_before_tone_onset;

%aline time to outcome
t_before_answer = 1000;
t_after_answer  = 4000;
Answertime_omiss = Answer_time(inds_use);
MinAnswerT = min(Answertime_omiss);
if MinAnswerT <= t_before_answer
    t_before_answer = MinAnswerT;
end
MinAnswer2endT = round(min(trial_all_time(inds_use)' - double(Answertime_omiss)));
if MinAnswer2endT <= t_after_answer
    t_after_answer = MinAnswer2endT;
end

for triN = 1:length(Answertime_omiss)
    Time_alined_outcome(triN,:) = double(Answertime_omiss(triN)) + [-t_before_answer+1 t_after_answer];
end
aline_time_window_outcome =  t_after_answer + t_before_answer;

end

