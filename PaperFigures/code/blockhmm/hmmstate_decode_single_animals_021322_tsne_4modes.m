%% classify and visualize state distributions
paths = pathsetup('matchingsim');

% Load classification info file
version = '121821bK4';
opts.rootdir = fullfile(paths.blockhmmfitpath, version);
opts.plotting = 1;
opts.save = 0;

classification_info_files = dir(fullfile(opts.rootdir, '*classification_info*_v1.mat'));
assert(numel(classification_info_files) == 1);
load(fullfile(classification_info_files(1).folder, classification_info_files(1).name),...
    'blockhmm_idx', 'statesFlat', 'folders', 'aggparams_native');

animalinfo = struct;

nstates_hmm = max(blockhmm_idx);
nstates_regimes = max(statesFlat);
assert(nstates_regimes == 6);
fprintf('Number of block hmm states = %d\n', nstates_hmm);


%% Perform the necessary state permutations for each animal
% We will permute the id's stored in zstates
% according to the id's in statesFlat
counter = 1; %keeps track of current position in statesFlat
animalinfo = struct;

for i = 1:numel(folders)
    load(fullfile(folders(i).folder, folders(i).name), 'zstates', 'params', 'lengths', 'transmat');
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
    
    [identityPermuted, stateSortIdx] = sort(statesFlat_extracted);
    transmatPermuted = transmat(stateSortIdx,stateSortIdx);
    
    labels = {};
    for j = 1:numel(identityPermuted)
        if identityPermuted(j) <= 4
            labels{j} = ['Q', num2str(identityPermuted(j))];    
        else
            labels{j} = ['IB', num2str(identityPermuted(j))];
        end
    end
    
    figure;
    imagesc(transmatPermuted)
    colormap gray
    caxis([0 1])
%     colorbar;
    axis xy
    mymakeaxis('x_label', 'State i', 'y_label', 'State i + 1', 'xticks', 1:4,...
        'xticklabels', labels, 'yticks', 1:4, 'yticklabels', labels, 'xytitle', animal_name);
    
    counter = counter + n_zstates;    
    
    % Break down into individual sessions
    assert(sum(lengths) == numel(zclassified));
    zclassified_splits = mat2cell(zclassified, 1, lengths);
    
    
    animalinfo(i).animal = animal_name;
    animalinfo(i).zclassified = zclassified_splits;
    animalinfo(i).classes = sort(statesFlat_extracted);
    
    filename = fullfile(paths.figpath, 'hmmblockFigs/transmat_animals/',...
        sprintf('%s_transmat_tsne_021322.pdf', animal_name));
    if ~exist(filename, 'file')
        saveas(gcf, filename);
    end

end

%% Plot identity of HMM modes per animal
hmmidentities = [];
namelst = {};
for i = 1:numel(animalinfo)
    hmmidentities(end+1,:) = animalinfo(i).classes;
    namelst{i} = animalinfo(i).animal;
end

meanmode = mean(hmmidentities, 2);
[~,idx] = sort(meanmode);
namelst = namelst(idx);
hmmidentities = hmmidentities(idx,:);

idx2 = 22 - [1,2,3,4,5,7,9,11,17,19,6,8,12,13,14,15,18,10,16,20,21];
idx2 = idx2(end:-1:1);
hmmidentities = hmmidentities(idx2,:);
namelst = namelst(idx2);

cmap = brewermap(6, 'Set1');
cmap = cmap([2,1,5,6,4,3],:);

figure;
imagesc(hmmidentities);
colormap(cmap);
hold on
vline((0:6) + 0.5, 'k');
hline((0:numel(animalinfo)) + 0.5, 'k');

axis xy
mymakeaxis('x_label', 'HMM mode', 'y_label', 'Animal', ...
    'font_size', 22, 'xticks', 1:6, 'yticks', 1:numel(namelst), 'yticklabels', namelst)








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
    xlimmax = min(nsess, 37);

    if opts.plotting
        % Plot bar graph to show the composition
        figure('Position', [440,423,736,375], 'Name', animalinfo(id).animal);
        h = bar(composition,'stacked');
        cols = paperaesthetics;
        colors = cols.colors;
%         colors = brewermap(6, 'Set1');
%         orders = [2, 1, 5, 4, 3];
        for i = 1:6
            h(i).FaceColor = colors(i,:);
            h(i).ShowBaseLine = 'off';
        end
        xlim([0.5 xlimmax + 0.5])
        ylim([0 1])
        mymakeaxis('x_label', 'Session #', 'y_label', 'Fraction', 'xticks', 0:5:xlimmax)
%         legend(h(1:5), {'Q-learning 1', 'Q-learning 2', 'Q-learning 3', 'Inference-based 4', 'Inference-based 5'}, 'Position', [0.4,0.42,0.1,0.1], ...
%             'FontSize', 14);
        filename = fullfile(paths.figpath, 'hmmblockFigs/compositions_animals/',...
            sprintf('%s_training_evolution_tsne_011622.pdf', animalinfo(id).animal));
        if ~exist(filename, 'file')
            saveas(gcf, filename);
        end
    end
end
% close all
close(f)


%% Save if requested
savefilename = fullfile(opts.rootdir, sprintf('hmmblock_composition_info_%s_v3', version));

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

paperaesthetics;
% colors = brewermap(6, 'Set1');
% colors = colors([1,2,3,4,5,6],:);
% colors(4,:) = [0,0,0];

figure;
hold on
lines = [];
for i =1:6
    sline = pad_to_same_length(extracted_all, i);
    N = sum(~isnan(sline));
    if i == 4
        h = errorbar(1:size(sline, 2), nanmean(sline, 1), nanstd(sline, [], 1) ./ sqrt(N), ...
            'o-', 'Color', 'k', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
    else
        h = errorbar(1:size(sline, 2), nanmean(sline, 1), nanstd(sline, [], 1) ./ sqrt(N), ...
            'o-', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    end
    lines(i) = h;
    xlim([1, 40])
    ylim([0, 1])
end

mymakeaxis('x_label', 'Session', 'y_label', 'Fraction', 'xticks', 0:5:40, 'font_size', 22)
l = legend(lines, {'Q1', 'Q2', 'Q3', 'Q4', 'IB5', 'IB6'});
l.Title.String = 'Regime';
l.Title.FontSize = 12;
l.FontSize = 12;


%%
sline = pad_to_same_length(extracted_all, 1);
means1 = nanmean(sline, 1);

sline = pad_to_same_length(extracted_all, 4);
means4 = nanmean(sline, 1);

sline = pad_to_same_length(extracted_all, 5);
means5 = nanmean(sline, 1);

sline = pad_to_same_length(extracted_all, 6);
means6 = nanmean(sline, 1);


