function [auc] = cal_auc2(neuroactivity, type)

%cal roc
threshold = unique(neuroactivity);
numpoint = length(threshold) + 1;

numpos = sum(type == 1);
numneg = sum(type == 0);

roc = nan(numpoint, 2);
for i = 1:numpoint-1
    pre_pos_idx = neuroactivity >= threshold(i);
    roc(i,1) = sum(type(pre_pos_idx)==1)/numpos;%true positive rate
    roc(i,2) = sum(type(pre_pos_idx)==0)/numneg;%false positive rate
end
roc(numpoint,1) = 0;
roc(numpoint,2) = 0; 

% %plot roc
% plot(roc(:,2), roc(:,1));
% xlabel('False positive rate');
% ylabel('True positive rate');

%cal auc
auc = -1*sum(diff(roc(:,2)).*(roc(1:end-1,1)+roc(2:end,1))/2);

if mean(neuroactivity(type == 1))<mean(neuroactivity(type == 0))
    auc = 1-auc;
end

end