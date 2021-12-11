% 9.29.21 expfit params updated
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_092921.mat');

% 9.30.21 prob expfit params
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_prob_093021.mat');

% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_100821_withknnMdl.mat')

% Load decoding models and performance of models (on simulated data)
% Models are located at Mdls1, Mdls09, Mdls08 etc (100-0, 90-10 and 80-20
% environments respectively).
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_101421_withknnMdl.mat')


% Load blockHMM results for a single animal
% This file contains info on: zstates and observations for all blocks of
% trials that are of interest
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/f02_hmmblockfit_102121c.mat')



%% visualize trials in state
allmeans = [];
for i = 1:4
    figure;
    obsfilt = obs(zstates == i-1, :);
    allmeans(i,:) = mean(obsfilt, 1);
    imagesc(obsfilt)
end

%% Fit allmeans to logistic curve




%%
% Calculate efficiency
effs = [];
for i = 1:4
    obsfiltered = obs(zstates == i-1,:);
    disp(size(obsfiltered))
    effs(i) = sum(obsfiltered(:) == 1) / numel(obsfiltered) / 15*20;
end


%%

% features_flat = [0.95  0.01  -2.7  -0.7780;
%     1  0.0434  -0.5266  -0.7780];

Mdl = Mdls1{1};
lapseFlat = params(3,:)';
effFlat = effs';
offsetFlat = params(1,:)';
slopesFlat = -params(2,:)';

% Decode for all animal sessions
% lapseFlat = reshape(explapses_all, [], 1);
% effFlat = reshape(expeff_all, [], 1);
% offsetFlat = reshape(expoffsets_all, [], 1);
% slopesFlat = reshape(expslopes_all, [], 1);

% Note: no normalization since normalization already handled in svm Mdl
features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];

features_flat(4, 2) = 0.07;
statesFlat = Mdl.predict(features_flat);
disp(statesFlat');
% eliminate nan's (nan's decode to 1)
% statesFlat(isnan(lapseFlat) | isnan(effFlat) | isnan(offsetFlat) | isnan(slopesFlat)) = nan;

%% classify and visualize state distributions
rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/';
folders = dir(fullfile(rootdir, '*hmmblockfit_102121.mat'));
load(fullfile(rootdir, 'allanimals_hmmblockfitsummary_102321.mat'));

filenames = {folders.name};
animals = cellfun(@(x) extractName(x), filenames, 'UniformOutput', false);


aggmeans = {};
aggparams = {};

animalinfo = struct;


for i = 1:numel(folders)
    load(fullfile(rootdir, folders(i).name));
    classes = classarr(i,:);
    zcopy = zstates;
    
    %permute zstates
    for istate = 0:3
        zcopy(zstates == istate) = classes(istate + 1);
    end
    
    % break back into lengths
    zstatecell = {};
    ctr = 1;
    for isess = 1:numel(lengths)
        zstatecell{isess} = zcopy(ctr:ctr + lengths(isess) - 1);
        ctr = ctr + lengths(isess);
    end
    
    animalinfo(i).zstatecell = zstatecell;

end


%% for each animal, plot the fraction of states on a rolling basis
window = 5; %sessions, actual window will be this value + 1
colors = brewermap(4, 'Set1');
colors(4,:) = [0 0 0];
colors = colors([2, 1, 3, 4],:); 

fracs_all_animals = {}; 

for id = 1:numel(animals)
    zcell = animalinfo(id).zstatecell;
    fracs_all = [];
    figure;
    hold on
    for i = 1:numel(zcell) - window
        cells = zcell(i:i+window);
        statearr = cell2mat(cells);
        disp(statearr)
        fracs = [];
        for istate = 1:4
            fracs(istate) = sum(statearr == istate) / numel(statearr);
        end
        fracs_all(i,:) = fracs;

    end
%     figure;
%     plot(fracs_all)

    lines = [];
    for i = 1:4
        l = plot(fracs_all(:,i), 'Color', colors(i,:));  
        lines(i) = l;
    end
    
    mymakeaxis('x_label', 'Sessions', 'y_label', 'Fraction of blocks', ...
        'xytitle', animals{id});
    lh = legend(lines, {'1', '2', '3', '4'}, 'Color', 'none');
    lh.Title.String = 'HMM mode';
    lh.Title.FontSize = 12;
    filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/hmmblockFigs/%s_trainingevolution.pdf',...
        animals{id});
%     saveas(gcf, filename);
    
    fracs_all_animals{id} = fracs_all;
    
%     title(animals{id})
end

% save('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/hmmblockfractions_102421b.mat', 'fracs_all_animals', 'animals');

%% Parse and average hmm summary fractions
rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/';
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/hmmblockfractions_102421b.mat');

% for i = 1:numel(animals)
%     filename = sprintf('%s_hmmblockfit_102121.mat', animals{i});
%     load(fullfile(rootdir, 'expdata', filename));
%     
%     assert(numel(fitrange) == size(fracs_all_animals{i}, 1) + 5);
% %     fprintf('%s: %d, %d\n', animals{i}, numel(fitrange), size(fracs_all_animals{i}, 1));
% end

animalranges = [6,29;
    6,29;
    7,28;
    6,29;
    7,20;
    6,20;
    6,22;
    8,22;
    2,16;
    8,27;
    5,48;
    3,26;
    6,20;
    6,15;
    7,23];
animalnames = {'f01', 'f02', 'f11', 'f12', 'f16', 'f17',...
    'f20', 'f21', 'fh02', 'fh03', 'e35', 'e54', 'e57', 'e46', 'e56'};
    
extracted_all = {};
for i = 1:numel(animalnames)
    pos = find(strcmp(animals, animalnames{i}));
    range = animalranges(i,:);
     
    
    filename = sprintf('%s_hmmblockfit_102121.mat', animalnames{i});
    load(fullfile(rootdir, 'expdata', filename));
    
    trialsrange = (max(range(1), fitrange(1))):(range(2)-5);
    
    trialidx = trialsrange - double(fitrange(1)) + 1;
    
    extracted = fracs_all_animals{pos}(trialidx,:);
    
    extracted_all{i} = extracted;
    
    
    fprintf('%s: fitrange: %d, %d, good range: %d, %d\n', animalnames{i}, fitrange(1), fitrange(end),...
        range(1), range(2));
%     disp(trialidx)
    
    
end

colors = brewermap(4, 'Set1');
colors = colors([2, 1, 3, 4],:);
colors(4,:) = [0,0,0];

% s1 = pad_to_same_length(extracted_all, 1);
% s2 = pad_to_same_length(extracted_all, 2);
% s3 = pad_to_same_length(extracted_all, 3);
% s4 = pad_to_same_length(extracted_all, 4);


figure;
hold on
lines = [];
for i =1:4
    sline = pad_to_same_length(extracted_all, i);
    h = errorbar(1:size(sline, 2), nanmean(sline, 1), nanstd(sline, [], 1) / sqrt(size(sline, 1)), ...
        'o-', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    lines(i) = h;
    xlim([1, 15])
    ylim([0, 1])
end

mymakeaxis('x_label', 'Session', 'y_label', 'Fraction', 'xticks', 0:5:15)
l = legend(lines, {'1', '2', '3', '4'});
l.Title.String = 'HMM mode';
l.Title.FontSize = 12;
l.FontSize = 12;


% errorbar(1:size(s2, 2), nanmean(s2, 1), nanstd(s2, [], 1)/ sqrt(size(s2, 1)));
% errorbar(1:size(s3, 2), nanmean(s3, 1), nanstd(s3, [], 1)/ sqrt(size(s3, 1)));
% errorbar(1:size(s4, 2), nanmean(s4, 1), nanstd(s4, [], 1)/ sqrt(size(s4, 1)));



function name = extractName(str)
splits = strsplit(str, '_');
name = splits{1};

end
















