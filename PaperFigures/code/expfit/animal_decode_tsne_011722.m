% Load expfit data
paths = pathsetup('matchingsim');
version = '122221b';
load(fullfile(paths.expdatapath, version, sprintf('expfit_params_%s.mat', version)));

% Load the knn model
opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10_tsne.mat';
opts.svmmodelpath = paths.svmmodelpath;
load(fullfile(opts.svmmodelpath, opts.model_name), 'Models');
modelid = 24;

Mdl = Models{modelid}{1};




%%
% Decode for all animal sessions
lapseFlat = reshape(explapses_all, [], 1);
effFlat = reshape(expeff_all, [], 1);
offsetFlat = reshape(expoffsets_all, [], 1);
slopesFlat = reshape(expslopes_all, [], 1);

% Note: no normalization since normalization already handled in svm Mdl
% features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];
features_flat = [-offsetFlat slopesFlat lapseFlat effFlat];

statesFlat = Mdl.predict(features_flat);
% eliminate nan's (nan's decode to 1)
statesFlat(isnan(lapseFlat) | isnan(effFlat) | isnan(offsetFlat) | isnan(slopesFlat)) = nan;


% Unflatten
states = reshape(statesFlat, size(expeff_all, 1), []);

%%
statecounts = [];
nstates = nanmax(states(:));
for i = 1:nstates
    statecounts(i,:) = sum(states == i, 1);
end

cols = paperaesthetics;
colors = cols.colors;

statefrac = statecounts ./ sum(statecounts);
statecumulative = cumsum(statefrac);

figure('Position', [440,423,736,375]);
h = bar(statefrac','stacked');
for i = 1:nstates
    h(i).FaceColor = colors(i,:);
    h(i).ShowBaseLine = 'off';
end
xlim([0.5,30.5])
mymakeaxis('x_label', 'Training days', 'y_label', 'Fraction of animals', 'xticks', 0:10:30,...
    'font_size', 22)
l = legend(h(1:nstates), {'Q1', 'Q2', 'Q3', 'Q4'}, 'Position', [0.4,0.42,0.1,0.1], ...
    'FontSize', 14);

%%
figure;
errorbar(1:size(states, 2), nanmean(states), nanstd(states, [], 1) / sqrt(size(states, 1)), 'o',...
    'MarkerFaceColor', cols.bluecol, 'MarkerEdgeColor', cols.bluecol, 'Color', cols.bluecol,...
    'LineWidth', 2)
xlim([0, 30])
ylim([1 3])
mymakeaxis('x_label', 'Training days', 'y_label', 'Mean decoded state', ...
    'yticks', 1:3, 'font_size', 22)
