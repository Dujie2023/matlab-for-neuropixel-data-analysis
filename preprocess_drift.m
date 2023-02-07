%% load data

np_data_pathway
bin_name = 'mouse7_20221013_RSC2L_g0_t0.imec0.ap.bin';
extract_waveformdata_dj(np_data_pathway, bin_name, rez)
%% plot drift for all units
%x axis:time
%y axis:good units
%value: current channel - mode(all channel)

