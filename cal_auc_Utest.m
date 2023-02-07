function [auc] = cal_auc_Utest(neuroactivity, type)
%Using the Mann–Whitney U test(also called the Mann–Whitney–Wilcoxon(MWW),
%Wilcoxon rank-sum test, or Wilcoxon–Mann–Whitney test) calculate AUC

x = neuroactivity(type == 1);
y = neuroactivity(type == 0);

n1 = length(x);
n2 = length(y);

[~,~,stats] = ranksum(x,y);
U = stats.ranksum-n1*(n1+1)/2;
auc = U/(n1*n2);

if mean(x)<mean(y)
    auc = 1-auc;
end

end


