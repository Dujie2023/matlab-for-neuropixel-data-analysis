function [auc] = cal_auc(neuroactivity, type)
%input: neuroactivity double matrix (1, ntrials)
%input: type int32 matrix(1, ntrials)
%output: auc double

%cal roc
[neuroactivity_sorted, sort_idx] = sort(neuroactivity);
type_sorted = type(sort_idx);

numpos = sum(type_sorted == 1);
numneg = sum(type_sorted == 0);
numpoint =length(neuroactivity_sorted);

roc = nan(numpoint, 2);
for i = 1:numpoint
    roc(i,1) = sum(type_sorted(1:i) == 1)/numpos;%true positive rate
    roc(i,2) = sum(type_sorted(1:i) == 0)/numneg;%false positive rate
end

% %plot roc
% plot(roc(:,2), roc(:,1));
% xlabel('False positive rate');
% ylabel('True positive rate');

%cal auc
auc = sum(diff(roc(:,2)).*(roc(1:end-1,1)+roc(2:end,1))/2);
end

