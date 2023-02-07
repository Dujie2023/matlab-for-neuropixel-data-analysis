%load('SessAreaIndexData.mat');
%load('np_extracted.mat', 'UsedClus_IDs');
function [cluster_Ab, usedAb, flidx_Ab] = extract_brainarea(SessAreaIndexStrc, UsedClus_IDs)
UsedAbIdx = find(SessAreaIndexStrc.UsedAbbreviations);
usedAb = cell(1,length(UsedAbIdx));

SessAllAb = fields(SessAreaIndexStrc);
cluster_Ab = nan(length(UsedClus_IDs),1);
flidx_Ab_IDs = [];

for i = 1:length(UsedAbIdx)
    usedAb{i} = SessAllAb{UsedAbIdx(i)};
    temp = getfield(SessAreaIndexStrc,SessAllAb{UsedAbIdx(i)}); 
    iAbidx = temp.MatchUnitRealIndex;
    flidx_Ab_IDs = union(flidx_Ab_IDs, iAbidx);
    for j = 1:length(iAbidx)
        cluster_Ab(UsedClus_IDs == iAbidx(j)) = i;
    end
end

[~, flidx_Ab, ~] = intersect(UsedClus_IDs,flidx_Ab_IDs);
end