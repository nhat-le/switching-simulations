filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata';
files = dir(fullfile(filedir, '*hmmblockfit_102121.mat'));

animals = {};
ranges = {};
for i = 1:numel(files)
    parts = strsplit(files(i).name, '_');
    animals{i} = parts{1};
    
    load(fullfile(files(i).folder, files(i).name));
    ranges{i} = fitrange;
end

save(fullfile(filedir, 'fitranges_102121.mat'), 'animals', 'ranges')