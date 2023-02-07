% read targetbrainregions.xlsx
% save corresponding index of abstract regions to a struct
txt = readcell("TargetbrainRegions.xlsx");
infotable = cell2table(txt(2:end,:));
infotable.Properties.VariableNames = txt(1,:);

abarea = [];
isabarea = true(size(infotable,1),1);
for i = 1:size(infotable,1)
    abidx = infotable{i,1};
    isabarea(i) = ismissing(string(abidx{:}));
end

abinfotable = infotable(~isabarea,:);
numabarea = size(abinfotable,1);
abinfo =  table2struct(abinfotable);
for i = 1:size(abinfotable,1)
    abidx_1 = strsplit(string(abinfotable{i,1}),';');
    abidx_temp = [];
    abidx_temp1 = [];
    for j = 1:size(abidx_1,2)
        [abidx_2,matches] = strsplit(abidx_1(j),'-');
        if ~ismissing(matches)
            abidx_temp1 = str2num(abidx_2(1)):1:str2num(abidx_2(2));
        else
            abidx_temp1 = str2num(abidx_2);
        end
        abidx_temp = [abidx_temp abidx_temp1];
    end
    abinfo(i).Index = int16(abidx_temp);
end
save('AbstractBrainRegions.mat', 'abinfo');