expfitdate = '040122';

load(fullfile('optodata/', expfitdate, 'opto_hmm_info.mat'));

sesscounts = nan(numel(animalinfo), 2, 2);
flatcounts = nan(numel(animalinfo), 4);


for animalID = 1:numel(animalinfo)
    out = helper.get_flags(animalinfo, animalID);
    for i = 1:2
        for j = 1:2
            counts = sum(out.power_flags{i} & ...
                out.period_flags{j} & sum(cell2mat(out.area_flags'), 1));
            sesscounts(animalID, i, j) =  counts;
            flatcounts(animalID, (i - 1) * 2 + j) = counts; 
        end
    end
end


%%
figure;
imagesc(flatcounts);
axis xy
hold on

for i = 1:4
    for j = 1:4
        if (j == 3 && i == 4) || (j == 4 && i < 4)
            color = 'k';
        else
            color = 'w';
        end
        
        text(j,i, num2str(flatcounts(i,j)), 'Color', color, 'FontSize', 20)
        
    end
end


mycolormap = customcolormap([0 1], {'#2166ac','#ffffff'});
xticks(1:4);
xticklabels({'LC', 'LO', 'HC', 'HO'})
yticks(1:4);
yticklabels({'f26', 'f27', 'f29', 'f32'})
colormap(mycolormap);

set(gca, 'FontSize', 20)












