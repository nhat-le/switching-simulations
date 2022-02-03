% Load prob expfit data: 9.30.21 prob expfit params (prob = 0.9, 0.8, 0.7, 6 animals)
paths = pathsetup('matchingsim');
load(fullfile(paths.expdatapath, 'expfit_params_prob_101421.mat'));

% Load the knn model
opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10b_tsne.mat';
opts.svmmodelpath = paths.svmmodelpath;
load(fullfile(opts.svmmodelpath, opts.model_name));
modelid = 24;

% Load prob expfit data
prob_animals = {'f11', 'f12', 'f16', 'f17', 'f20', 'f21'};
version = '122221b';
load(fullfile(paths.expdatapath, version, 'fitparams_session_averaged_122221b.mat'))

% For each animal, identify the prob flags
probflags_all = {};
explapses = {};
expoffsets = {};
expslopes = {};
expeffs = {};
expprobs = {};
for i = 1:numel(prob_animals)
    animal = prob_animals{i};
    load(fullfile(paths.expdatapath, version, sprintf('%s_all_sessions_%s.mat', ...
        animal, version)), 'probflags'); %'fitparams_session_averaged_122221b.mat'))
    probflags_all{end+1} = probflags;
    assert(numel(fitparams_all.(animal).eff) == numel(probflags));
    
    params = fitparams_all.(animal);
    [offsets, slopes, lapses, effs] = parse_params(params);
    
    explapses{i} = lapses(probflags > 0);
    expslopes{i} = slopes(probflags > 0);
    expoffsets{i} = offsets(probflags > 0);
    expeff{i} = effs(probflags > 0);
    expprobs{i} = probflags(probflags > 0); 
    
end


explapses_all = pad_to_same_length(explapses);
expoffsets_all = pad_to_same_length(expoffsets);
expslopes_all = pad_to_same_length(expslopes);
expeff_all = pad_to_same_length(expeff);
expprobs_all = pad_to_same_length(expprobs);


%%
% Decode for all animal sessions
lapseFlat = reshape(explapses_all, [], 1);
effFlat = reshape(expeff_all, [], 1);
offsetFlat = reshape(expoffsets_all, [], 1);
slopesFlat = reshape(expslopes_all, [], 1);

if exist('expprobs_all', 'var')
    probsFlat = reshape(expprobs_all, [], 1);
else
    probsFlat = ones(1, numel(slopesFlat));
end

% Note: no normalization since normalization already handled in svm Mdl
% features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];
features_flat = [-offsetFlat -slopesFlat lapseFlat effFlat];


probs = [1, 0.9 0.8 0.7];
mdls_all = {Mdl10{1}, Mdl09{1}, Mdl08{1}, Mdl07{1}};
statesFlat = nan(numel(probsFlat), 1);
for i = 1:numel(probs)
%     sessi = find(probsFlat == prob);
    featuresi = features_flat(probsFlat == probs(i), :);
    statesFlat(probsFlat == probs(i)) = mdls_all{i}.predict(featuresi);

end
    
% statesFlat = Mdl.predict(features_flat);
% eliminate nan's (nan's decode to 1)
statesFlat(isnan(lapseFlat) | isnan(effFlat) | isnan(offsetFlat) | isnan(slopesFlat)) = nan;


% Unflatten
states = reshape(statesFlat, size(expeff_all, 1), []);

%% Aggregate by probs and plot
states09 = statesFlat(probsFlat == 0.9);
states08 = statesFlat(probsFlat == 0.8);
states07 = statesFlat(probsFlat == 0.7);
nstates = nanmax(statesFlat);

for i = 1:nstates
    frac09(i) = sum(states09 == i) / numel(states09);
end

for i = 1:nstates
    frac08(i) = sum(states08 == i) / numel(states08);
end

for i = 1:nstates
    frac07(i) = sum(states07 == i) / numel(states07);
end

figure;
cols = paperaesthetics;
colors = cols.colors;

h = bar([frac09; frac08; frac07], 'stack');
for i = 1:nstates
    h(i).FaceColor = colors(i,:);
    h(i).ShowBaseLine = 'off';
end
xlim([0.5,4.5])
labels = {'Q1', 'Q2', 'Q3', 'Q4', 'IB5', 'IB6'};
mymakeaxis('x_label', 'Probability', 'y_label', 'Fraction of sessions', 'xticks', 1:3, ...
    'xticklabels', {'90-10', '80-20', '70-30'}, 'font_size', 22)
l = legend(h(1:nstates), labels(1:nstates), 'Position', [0.78,0.38,0.15,0.17], ...
    'FontSize', 14);



function [offsets, slopes, lapses, effs] = parse_params(params)

slopesL = params.pL(:,1);
slopesR = params.pR(:,1);
offsetsL = params.pL(:,2);
offsetsR = params.pR(:,2);
lapsesL = params.pL(:,3);
lapsesR = params.pR(:,3);

slopes = nanmin([slopesL, slopesR], [], 2)';
offsets = nanmean([offsetsL, offsetsR], 2)';
lapses = nanmax([lapsesL, lapsesR], [], 2)';
effs  = params.eff;




end
