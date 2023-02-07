function [sig_thre05, sig_thre01] = auc_sigthreshold(all_neuroactivity, type, repeat_times)
%input
%all_neuroactivity: neuroactivity of all neurons
%type: label, e.g. blocktype, tritype

%output
%sigthre05: 0.05
%sigthre01: 0.01

numneuron = size(all_neuroactivity, 1);
AUC_random = zeros(repeat_times, numneuron);
numlabel = length(type);
for i = 1:repeat_times
    type_i = type(randperm(numlabel));
    fprintf(['repeat time : ', num2str(i), '\n']);
    for jneuron = 1:numneuron
        AUC_random(i, jneuron) = cal_auc2(all_neuroactivity(jneuron,:), type_i);
    end
end

%distribution
AUC_random = AUC_random(:);
figure;
histogram(AUC_random);

sig_thre05 = prctile(AUC_random,95);
sig_thre01 = prctile(AUC_random,99);
end

