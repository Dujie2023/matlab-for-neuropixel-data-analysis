function [sig_thre05, sig_thre01] = parfor_auc_sigthreshold(all_neuroactivity, type, repeat_times, savefolder, varname,useUtest)
%input
%all_neuroactivity: neuroactivity of all neurons
%type: label, e.g. blocktype, tritype

%output
%sigthre05: 0.05
%sigthre01: 0.01

CoreNum=6; %设定机器CPU核心数量
if isempty(gcp('nocreate')) %如果并行未开启
    parpool(CoreNum);
end

numneuron = size(all_neuroactivity, 1);
AUC_random = zeros(repeat_times, numneuron);
numlabel = length(type);

fprintf(['calculating AUC significant threshold... ', '\n']);
if useUtest == 1
    parfor i = 1:repeat_times
        type_i = type(randperm(numlabel));
        %fprintf(['repeat time : ', num2str(i), '\n']);
        for jneuron = 1:numneuron
            AUC_random(i, jneuron) = cal_auc_Utest(all_neuroactivity(jneuron,:), type_i);
        end
    end
else
    parfor i = 1:repeat_times
        type_i = type(randperm(numlabel));
        %fprintf(['repeat time : ', num2str(i), '\n']);
        for jneuron = 1:numneuron
            AUC_random(i, jneuron) = cal_auc2(all_neuroactivity(jneuron,:), type_i);
        end
    end
end

%distribution
AUC_random = AUC_random(:);
figure;
histogram(AUC_random);
xlabel('AUC');
saveas(gcf, fullfile(savefolder,[varname 'AUC sig histogram.png' ]));
saveas(gcf, fullfile(savefolder,[varname 'AUC sig histogram.fig' ]));
close;

sig_thre05 = prctile(AUC_random,95);
sig_thre01 = prctile(AUC_random,99);

end

