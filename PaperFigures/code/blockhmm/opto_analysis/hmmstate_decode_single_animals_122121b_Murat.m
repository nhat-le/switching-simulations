%% classify and visualize state distributions
paths = pathsetup('matchingsim');

% Load classification info file
version = '012122_Murat';
opts.rootdir = fullfile(paths.blockhmmfitpath, version);
opts.plotting = 1;
opts.save = 1;

classification_info_files = dir(fullfile(opts.rootdir, '*classification_info*.mat'));
assert(numel(classification_info_files) == 1);
load(fullfile(classification_info_files(1).folder, classification_info_files(1).name),...
    'blockhmm_idx', 'statesFlat', 'folders', 'aggparams_native');

animalinfo = struct;

nstates_hmm = max(blockhmm_idx);
nstates_regimes = max(statesFlat);
assert(nstates_regimes == 5);
fprintf('Number of block hmm states = %d\n', nstates_hmm);


%% Perform the necessary state permutations for each animal
% We will permute the id's stored in zstates
% according to the id's in statesFlat
counter = 1; %keeps track of current position in statesFlat
animalinfo = struct;

for i = 1:numel(folders)
    load(fullfile(folders(i).folder, folders(i).name), 'zstates', 'params', 'lengths');
    parts = strsplit(folders(i).name, '_');
    animal_name = parts{1};
    zclassified = zstates;
    
    n_zstates = max(zstates) + 1;
    assert(n_zstates == size(params, 2))
    
    % make sure the states are in the right order
    assert(sum(sum(aggparams_native{i}(1:3,:) ~= params)) == 0);
    
    statesFlat_extracted = statesFlat(counter : counter + n_zstates - 1);
    for istate = 1:n_zstates
        zclassified(zstates == istate - 1) = statesFlat_extracted(istate);
    end
    
    
    counter = counter + n_zstates;    
    
    % Break down into individual sessions
    assert(sum(lengths) == numel(zclassified));
    zclassified_splits = mat2cell(zclassified, 1, lengths);
    
    
    animalinfo(i).animal = animal_name;
    animalinfo(i).zclassified = zclassified_splits;

end

%% Plot state evolution profile for each animal
f = waitbar(0);
for id = 1:numel(animalinfo)
    waitbar(id / numel(animalinfo), f);
    zclassified = animalinfo(id).zclassified;
    composition = nan(numel(zclassified), nstates_regimes);
    for i = 1:numel(zclassified)
        z_single_session = zclassified{i};
        for j = 1:nstates_regimes
            composition(i,j) = sum(z_single_session == j) / numel(z_single_session);   
        end
    end
    
    animalinfo(id).composition = composition;
    nsess = size(composition, 1);
    xlimmax = min(nsess, nsess);

    if opts.plotting
        % Plot bar graph to show the composition
        figure('Position', [440,423,736,375], 'Name', animalinfo(id).animal);
        h = bar(composition,'stacked');
        colors = brewermap(6, 'Set1');
        orders = [2, 1, 5, 4, 3];
        for i = 1:5
            h(i).FaceColor = colors(orders(i),:);
            h(i).ShowBaseLine = 'off';
        end
        xlim([0.5 xlimmax + 0.5])
        mymakeaxis('x_label', 'Session #', 'y_label', 'Fraction', 'xticks', 0:5:30)
%         legend(h(1:5), {'Q-learning 1', 'Q-learning 2', 'Q-learning 3', 'Inference-based 4', 'Inference-based 5'}, 'Position', [0.4,0.42,0.1,0.1], ...
%             'FontSize', 14);
        filename = fullfile(paths.figpath, 'hmmblockFigs/compositions_animals/',...
            sprintf('%s_training_evolution.pdf', animalinfo(id).animal));
%         saveas(gcf, filename);
    end
end
% close all
close(f)


%% Save if requested
savefilename = fullfile(opts.rootdir, sprintf('hmmblock_composition_info_%s', version));

if opts.save && ~exist(savefilename, 'file')
    save(savefilename, 'opts', 'animalinfo');
    fprintf('File saved!\n');
else
    fprintf('Skipping save...\n');   
end


%% Parse and average hmm summary fractions
extracted_all = {};

for i = 1:numel(animalinfo)
    animalname = animalinfo(i).animal;
    filename = sprintf('%s_hmmblockfit_*%s.mat', animalname, version);
    blockfitfile = dir(fullfile(paths.blockhmmfitpath, version, filename));
    assert(numel(blockfitfile) == 1);
    load(fullfile(blockfitfile(1).folder, blockfitfile(1).name), 'fitrange');
    
    
    trialidx = 1:min(40, size(animalinfo(i).composition, 1));
    
    extracted = animalinfo(i).composition(trialidx,:);
    
    extracted_all{i} = extracted;
    
    
    
end

colors = brewermap(6, 'Set1');
colors = colors([2, 1, 5, 4, 3],:);
% colors(4,:) = [0,0,0];

figure;
hold on
lines = [];
for i =1:5
    sline = pad_to_same_length(extracted_all, i);
    N = sum(~isnan(sline));
    h = errorbar(1:size(sline, 2), nanmean(sline, 1), nanstd(sline, [], 1) ./ sqrt(N), ...
        'o-', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    lines(i) = h;
    xlim([1, 30])
    ylim([0, 1])
end

mymakeaxis('x_label', 'Session', 'y_label', 'Fraction', 'xticks', 0:5:30)
l = legend(lines, {'Q1', 'Q2', 'Q3', 'IB4', 'IB5'});
l.Title.String = 'Regime';
l.Title.FontSize = 12;
l.FontSize = 12;





