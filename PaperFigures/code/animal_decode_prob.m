% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/svm_Mdl.mat');
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params.mat');

% Decoding for probabilistic environments

% 9.29.21 expfit params updated
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_092921.mat');

% 9.30.21 prob expfit params
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_prob_093021.mat');

load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_092921_withMdl.mat')

probdecode = 1;

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
features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];

probs = [1, 0.9 0.8 0.7];
mdls_all = {Mdls1{1} Mdls09{1}, Mdls08{1}, Mdls07{1}};
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


%%
paperaesthetics;
errorbar(1:size(states, 2), nanmean(states), nanstd(states, [], 1) / sqrt(size(states, 1)), 'o',...
    'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol, 'Color', bluecol,...
    'LineWidth', 2)
xlim([0, 20])
ylim([1 4])
mymakeaxis('x_label', 'Training days', 'y_label', 'Mean decoded state')


%% Visualize animal behavior in mds space 
[res1, ~] = load_and_run(0);
[res2, ~] = load_and_run(0.1);
[res3, ~] = load_and_run(0.2);
[res4, ~] = load_and_run(0.3);
res_all = {res1, res2, res3, res4};


features_flatcopy = features_flat;
features_flatcopy(:,3) = -features_flat(:,3);
features_flatcopy(features_flatcopy < -20) = 3;

Yall = {};
for i = 1:numel(probs)
    featuresi = features_flatcopy(probsFlat == probs(i), :);
    features_norm = (featuresi - nanmean(res_all{i}.features, 1)) ./ nanstd(res_all{i}.features, [], 1);
    features_proj = features_norm * res_all{i}.V;
    Yall{i} = features_proj(:,1:2);
    disp(size(Yall{i}));
end


%% plot
figure;
hold on
for i = 1:numel(probs)
    subplot(2,2,i)
    hold on
    for j = 1:numel(unique(res_all{i}.idx))
        plot(res_all{i}.Y(res_all{i}.idx == j, 1), res_all{i}.Y(res_all{i}.idx == j, 2), '.')
    end

%     plot(Yall{i}(i,1,1), Yall{i}(i,1,2), 'ro')
    plot(Yall{i}(:, 1), Yall{i}(:,2), 'rx')
end
% ylim([-2 2])



%% 
for i = 1:15
    subplot(3,5,i)
    plot(states(i,:), 'o')
    ylim([0, 5])
    
end