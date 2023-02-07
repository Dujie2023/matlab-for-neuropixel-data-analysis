function waveform_feature_funcollect = waveform_funcollect
%function collections for calculating waveform features

waveform_feature_funcollect.cal_waveform_feature = @cal_waveform_feature;
waveform_feature_funcollect.cal_waveform_spread  = @cal_waveform_spread;
waveform_feature_funcollect.cal_wave_amplitude = @cal_wave_amplitude;
waveform_feature_funcollect.cal_wave_snr = @cal_wave_snr;

end

%% function for extracting waveform feature for a single cluster
function [wave_amplitude, trough2peakT, repolarT, prepeak2pospeakratio, troughpeakidx] = ...
                                cal_waveform_feature(AvgWaves,WaveWinSamples)
%feature 1 wave_amplitute: the difference (in microvolts) between the peak and 
%                          trough of the waveform on a single
%                          channel;according to Joshua H. Siegle, 2021
%feature 2 trough2peakT: the distance between the global minimum of the 
%                       curve and the following local maximum. 
%feature 3 repolarT: the distance between the late positive peak and the
%                    inflection point of the falling branch of the curve;
%                    according to Trainito et al., 2019
%feature 4: prepeak2postpeakratio

wave_amplitude = nan;
trough2peakT = nan;
repolarT = nan;
prepeak2pospeakratio = nan; 
troughpeakidx = [nan, nan];

%trough defined as minium point
[trough_am, trough_idx] = min(AvgWaves,[],'omitnan');
trough_idx = trough_idx(1);
trough_am = trough_am(1);

af_trough_Wave = AvgWaves(trough_idx:end);
[peak_am, peak_idx] = max(af_trough_Wave);
peak_am = peak_am(1);
peak_idx = peak_idx(1);

trough2peakT = peak_idx - 1;%/30 ms
wave_amplitude = peak_am - trough_am;
post_peak_idx = peak_idx + trough_idx - 1;

af_peak_Wave = AvgWaves(post_peak_idx:end);
if length(af_peak_Wave) >= 3    
    %find the inflection point
    d1 = gradient(af_peak_Wave);
    [~, inflepoint_idx] = min(d1);
    repolarT = inflepoint_idx(1) - 1;
end
    
%whether peaks before trough exist
try
    [bf_peak, ~,~,p] = findpeaks(AvgWaves(1:trough_idx));
    if ~isempty(bf_peak) 
        try
            [~, af_peaks_idx, ~, pafs] = findpeaks(af_trough_Wave);
            temp = find(af_peaks_idx == peak_idx);
            paf = pafs(temp(1));%for calculating pre2postpeakratio
            prepeak2pospeakratio = max(p)/paf(1);
        end
    else
        prepeak2pospeakratio = 0;
    end
end

troughpeakidx = [trough_idx, post_peak_idx];

end

%% function for calculating waveform spread for a single cluster
function [wave_spread] = cal_waveform_spread(AvgAllChnWaves, cChannel)
%wave_spread: Spatial extent (in ¦Ìm) of channels in which the waveform amplitude
%             exceeds 12% of the peak amplitude;;according to Joshua H. Siegle, 2021
%output wave_spread: unit channel_numbel
%AvgAllChnWaves:matrix(numChannel,diff(WaveWinSamples))

cl_range = [cChannel-100; cChannel+100];
if cl_range(1)<1
    cl_range(1) = 1;
end
if cl_range(2)>384
    cl_range(2) = 384;
end

am = zeros(diff(cl_range)+1,1);
for i = cl_range(1):cl_range(2)
    am(i-cl_range(1)+1) = cal_wave_amplitude(AvgAllChnWaves(i,:));
end
norm_ampltitude = am/max(am);

wave_spread = numel(find(norm_ampltitude>=0.12));

end
%% function for calculating waveform amplitude
function [ wave_amplitude] = cal_wave_amplitude(AvgWaves)
%trough defined as minium point
[trough_am, trough_idx] = min(AvgWaves,[],'omitnan');
trough_idx = trough_idx(1);
trough_am = trough_am(1);

[af_peak, ~] = max(AvgWaves(trough_idx:end));
wave_amplitude = af_peak(1) - trough_am;

end

%% function for calculating signal-to-noise ratio
% defined as the ratio between the waveform amplitude 
%   and 2¡Á the standard deviation of the residual waveforms73

function [wave_SNR] = cal_wave_snr(channel_waveform)
avg_wave = mean(channel_waveform, 'omitnan');
am = max(avg_wave)-min(avg_wave);
e = channel_waveform - repmat(avg_wave, size(channel_waveform, 1), 1);
wave_SNR = am/(2*std(e(:),'omitnan'));
end