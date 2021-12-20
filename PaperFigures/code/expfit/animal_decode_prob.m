% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/svm_Mdl.mat');
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params.mat');

% Decoding for probabilistic environments

% 9.29.21 expfit params updated (prob = 1, 15 animals)
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_092921.mat');

% 9.30.21 prob expfit params (prob = 0.9, 0.8, 0.7, 6 animals)
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_prob_101421.mat');

load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_100821_withknnMdl.mat')

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

usepca = 0;
[res1, ~] = load_and_run(0, 0);
[res2, opts] = load_and_run(0.1, 0);
[res3, ~] = load_and_run(0.2, 0);
[res4, ~] = load_and_run(0.3, 0);

% [res1pca, ~] = load_and_run(0, 1);
% [res2pca, ~] = load_and_run(0.1, 1);
% [res3pca, ~] = load_and_run(0.2, 1);
% [res4pca, ~] = load_and_run(0.3, 1);


% Split/combine so that we end up with 5 clusters for each space
% For prob = 1, will split cluster 4 into two, based on results of clusters
% for prob = 0.8
idxprob1 = res1.idx;
idxprob08 = res3.idx;

idxprob1(idxprob1 == 4 & idxprob08 == 5) = 5;
% [idxQ, idxIB] = reshapeidx(idxprob1, res1);

% for res1, split cluster 4 into two, as done above
res1new = res1;
res1new.idx = idxprob1;

% for res2, combine 4 and 5
res2new = res2;
idxprob09 = res2.idx;
idxprob09(idxprob09 == 5) = 4;
idxprob09(idxprob09 == 6) = 5;
res2new.idx = idxprob09;

res3new = res3;

res4new = res4;

features_proj = {};
res_all = {res1new, res2new, res3new, res4new};


%%
i = 3;
%
load(sprintf('%s/svmresults_from_pickle_092221_prob%.2f.mat', opts.rootdir, 1 - probs(i)));

[idxQ, idxIB] = reshapeidx(res1new.idx, res1new);


idxQrep = repmat(idxQ, [1 1 50]);
idxIBrep = repmat(idxIB, [1 1 50]);

%unroll the matrices
idxQall = reshape(idxQrep, [], 1);
idxIBall = reshape(idxIBrep, [], 1);

IBeffall = reshape(IBeff_arr, [], 1);
IBslopeall = reshape(IBslope_arr, [], 1);
IBlapseall = reshape(IBlapse_arr, [], 1);
IBoffsetall = reshape(IBoffset_arr, [], 1);
Qeffall = reshape(Qeff_arr, [], 1);
Qslopeall = reshape(Qslope_arr, [], 1);
Qlapseall = reshape(Qlapse_arr, [], 1);
Qoffsetall = reshape(Qoffset_arr, [], 1);

Qoffsetall(Qoffsetall < -20) = 3;
IBoffsetall(IBoffsetall < -20) = 3;

features = [IBeffall IBlapseall -IBslopeall IBoffsetall;
    Qeffall Qlapseall -Qslopeall Qoffsetall];

% features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

features_norm = (features - mean(res_all{i}.features, 1)) ./ (std(res_all{i}.features, [], 1));
features_proj = features_norm * res_all{i}.V;



figure;
hold on
for j = 1:5
    plot(features_proj(labels == j, 1), features_proj(labels == j, 2), '.');
end
xlim([-5, 5])
ylim([-5, 5])

features_norm2 = (res_all{i}.features - mean(res_all{i}.features, 1)) ./ (std(res_all{i}.features, [], 1));
features_proj2 = features_norm2 * res_all{i}.V;
labels2 = res_all{i}.idx;
figure;
hold on
for j = 1:5
    plot(features_proj2(labels2 == j, 1), features_proj2(labels2 == j, 2), '.');
end


xlim([-5, 5])
ylim([-5, 5])

%%
figure;
handles = [];
for j = 1:4
    h1 = subplot(4,2,2*j - 1)
    plot(res_all{i}.features(res_all{i}.idx == 1, j));

    h2 = subplot(4,2,2*j)
    plot(features(labels == 1, j));
    linkaxes([h1 h2], 'y')
    handles = [handles h1 h2];
end
linkaxes(handles(1:2:end), 'x')
linkaxes(handles(2:2:end), 'x')



%%

features_proj = {};
for i = 1:numel(probs)
%     sessi = find(probsFlat == prob);
    featuresi = features_flat(probsFlat == probs(i), :);
    featuresi(:,3) = -featuresi(:,3);
    disp(size(featuresi))
%     statesFlat(probsFlat == probs(i)) = mdls_all{i}.predict(featuresi);
    features_norm = (featuresi - nanmean(res_all{i}.features, 1)) ./ nanstd(res1.features, [], 1);
    features_proj{i} = features_norm * res_all{i}.V;
end

% labels = res1.idx;

% features_norm = (features_flat - nanmean(res1.features, 1)) ./ nanstd(res1.features, [], 1);
% features_proj = features_norm * res1.V;
% Yanimals = reshape(features_proj(:, 1:2), size(expeff_all, 1), [], 2);



%% Aggregate by probs and plot
states09 = statesFlat(probsFlat == 0.9);
states08 = statesFlat(probsFlat == 0.8);
states07 = statesFlat(probsFlat == 0.7);

for i = 1:5
    frac09(i) = sum(states09 == i) / numel(states09);
end

for i = 1:5
    frac08(i) = sum(states08 == i) / numel(states08);
end

for i = 1:5
    frac07(i) = sum(states07 == i) / numel(states07);
end

figure;
% focus on f12 for main figure
colors = brewermap(6, 'Set1');
colorsbold = brewermap(6, 'Set1');

colors = colors([2,1,5,4,3],:);
colorsbold = colorsbold([2,1,5,4,3],:);
h = bar([frac09; frac08; frac07], 'stack');
for i = 1:5
    h(i).FaceColor = colors(i,:);
    h(i).ShowBaseLine = 'off';
end
xlim([0.5,4.5])
mymakeaxis('x_label', 'Probability', 'y_label', 'Fraction of sessions', 'xticks', 1:3, ...
    'xticklabels', {'90-10', '80-20', '70-30'})
l = legend(h(1:3), {'Class 1', 'Class 2', 'Class 3'}, 'Position', [0.78,0.38,0.15,0.17], ...
    'FontSize', 14);

%% Plot mds
clustid = 4; 
xmin = -6;
xmax = 4;
ymin = -4;
ymax = 4;

colorsbold = brewermap(6, 'Set1');
colorsbold = colorsbold([2,1,5,4,3],:);
colors = brewermap(6, 'Pastel1');
colors = colors([2,1,5,4,3],:);
figure;
hold on

simfeatures = res_all{clustid}.features;
Vmat = res_all{clustid}.V;
simfeatures_norm = (simfeatures - nanmean(simfeatures, 1)) ./ nanstd(simfeatures, [], 1);
simfeatures_proj = simfeatures_norm * Vmat;

xcoords = linspace(xmin, xmax, 100);
ycoords = linspace(ymin, ymax, 100);
% [xx, yy] = meshgrid(xcoords, ycoords, ycoords);

domains = zeros(numel(xcoords), numel(ycoords));

% for each (x,y), find the Y entry that is closest
for i = 1:numel(xcoords)
    for j = 1:numel(ycoords)
        xval = xcoords(i);
        yval = ycoords(j);
        
        D = [xval yval] - simfeatures_proj(:,1:2);
        
        dist = D(:,1).^2 + D(:,2).^2;
        
%         Mdls1{1}.predict([xval yval]
        
        id = argmin(dist);
        
        domains(i,j) = res_all{clustid}.idx(id);
        
        
    end
end


imagesc(domains', 'XData', xcoords, 'YData', ycoords);
for i = 1:5
    plot(simfeatures_proj(res_all{clustid}.idx == i, 1), simfeatures_proj(res_all{clustid}.idx == i, 2), 'o', ...
      'MarkerFaceColor', colorsbold(i,:), 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
  
end

datafeatures = features_flat(probsFlat == probs(clustid), :);
datafeatures(:,3) = -datafeatures(:,3);
datafeatures_norm = (datafeatures  - nanmean(simfeatures, 1)) ./ nanstd(simfeatures, [], 1);
datafeatures_proj = datafeatures_norm * Vmat;

% plot(Yanimals(1,1), Yanimals(1,2), 'rx', 'LineWidth', 3, 'MarkerSize', 15)
states_all = {states09, states08, states07};
statesarr = states_all{clustid - 1};

for i = 1:5
    plot(datafeatures_proj(statesarr == i, 1), datafeatures_proj(statesarr == i, 2), 'x', ...
      'MarkerFaceColor', colorsbold(i,:), 'MarkerEdgeColor', colorsbold(i,:), 'MarkerSize', 7,...
      'LineWidth', 2);
  
end


% plot(datafeatures_proj(statesarr == 1, 1), datafeatures_proj(statesarr == 1,2), 'bx', 'MarkerSize', 15, 'LineWidth', 1)
% plot(datafeatures_proj(statesarr == 2, 1), datafeatures_proj(statesarr == 2,2), 'rx', 'MarkerSize', 15, 'LineWidth', 1)
% plot(datafeatures_proj(statesarr == 3, 1), datafeatures_proj(statesarr == 3,2), 'yx', 'MarkerSize', 15, 'LineWidth', 1)
% plot(datafeatures_proj(statesarr == 4, 1), datafeatures_proj(statesarr == 4,2), 'kx', 'MarkerSize', 15, 'LineWidth', 1)
% plot(datafeatures_proj(statesarr == 5, 1), datafeatures_proj(statesarr == 5,2), 'gx', 'MarkerSize', 15, 'LineWidth', 1)


colormap(colors)
% title(all_animals(i,:))

xlim([xmin, xmax])
ylim([ymin, ymax])
mymakeaxis('x_label', 'PC1', 'y_label', 'PC2')


%% Visualize data features
figure;
handles = [];
for j = 1:4
    h1 = subplot(4,2,2*j - 1);
    plot(datafeatures(:, j));

    h2 = subplot(4,2,2*j);
    plot(simfeatures(:, j));
    linkaxes([h1 h2], 'y')
    handles = [handles h1 h2];
end
linkaxes(handles(1:2:end), 'x')
linkaxes(handles(2:2:end), 'x')


%%
% classifying point 9
dpoint = datafeatures_norm(9,:);
c1 = simfeatures_norm(167,:);
c2 = simfeatures_norm(17,:);

apoint = datafeatures(9,:);
a1 = simfeatures(167,:);
a2 = simfeatures(17,:);

rpoint = datafeatures_proj(9,:);
r1 = simfeatures_proj(167,:);
r2 = simfeatures_proj(17,:);



d1 = sum((dpoint - c1).^2);
d2 = sum((dpoint - c2).^2);


da1 = sum((apoint - a1).^2);
da2 = sum((apoint - a2).^2);

dr1 = sum((rpoint - r1).^2);
dr2 = sum((rpoint - r2).^2);









