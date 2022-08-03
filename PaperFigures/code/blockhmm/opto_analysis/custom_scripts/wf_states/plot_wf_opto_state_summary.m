folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/optodata/061522defaultK6';
filename = 'wf_hmm_info.mat';
load(fullfile(folder, filename))

behaviorfolder = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/processed/raw/behavior';
dirlst = dir(behaviorfolder);

animal_lst = {};
states_lst = {};

for i = 4:numel(dirlst)
    dirname = dirlst(i).name;
    parts = strsplit(dirname, '_');
    animal = parts{end};
    datestring = parts{1};

    animal_idx = strcmpi({animalinfo.animal}, animal);
    assert(sum(animal_idx) == 1);

    sessids = animalinfo(animal_idx).sessnames;
    filename = sprintf('%s_1_%s_block.mat', datestring, animal);
    found = 0;
    found_idx = -1;
    for sessid = 1:size(sessids, 1)
        if strcmpi(sessids(sessid,:), filename)
%             disp('found')
            found_idx = sessid;
            found = 1;
        end
    end
    assert(found)

    zstates = animalinfo(animal_idx).zclassified{found_idx};
    animal_lst{end+1} = animal;
    states_lst{end+1} = zstates;
    
end

%% Summarize state information per animal
animals = unique(animal_lst);
compositions = [];
for i = 1:numel(animals)
    states_all = states_lst(strcmp(animal_lst, animals{i}));
    states_all = cell2mat(states_all);
    for j = 1:6
        compositions(i,j) = sum(states_all == j);
    end
end

compositions = compositions ./ sum(compositions, 2);

figure('Position', [440,409,748,389]);
colors = brewermap(6, 'Set1');
colors = colors([2,1,5,6,4,3],:);
xvals = 1:6;
h = bar(xvals, compositions, 'stacked');

for i = 1:6
    h(i).FaceColor = colors(i,:);

end

mymakeaxis('x_label', '', 'y_label', 'Fraction of trials', ...
    'xticks', xvals, 'xticklabels', animals)













