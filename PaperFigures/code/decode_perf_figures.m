% Figures for decoding performance
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/svmresults_from_pickle_083121.mat');

mindecodingQ = min(decoding_perf, [], [3, 4]);
mindecodingIB = squeeze(min(decoding_perf, [], [1, 2]));

%%
produce_heatmap(mindecodingQ', epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Decoding performance', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');


%%
produce_heatmap(mindecodingIB', prewlst, pswitchlst, 'clim', [0.5, 1], 'legendname', 'Decoding performance', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);


%% Regime demarcation (with k-means clustering)
% Form the feature vectors
filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat';
load(filedir);

% rng(19)
% rng(20)
% rng(26)
% rng(29)
% rng(31)
% rng(123) %when N =5
rng(123)
N = 5;
Qeff_flat = reshape(efflist, [], 1);
Qlapse_flat = reshape(LapseL, [], 1);
Qslope_flat = reshape(PLslopelist, [], 1);
Qoffset_flat = reshape(PLoffsetlist, [], 1);

[Qxdim, Qydim] = size(efflist);

filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-withCorr-prob0.00to1.00-072121.mat';
load(filedir);

IBeff_flat = reshape(efflist, [], 1);
IBlapse_flat = reshape(LapseL, [], 1);
IBslope_flat = reshape(PLslopelist, [], 1);
IBoffset_flat = reshape(PLoffsetlist, [], 1);


[IBxdim, IBydim] = size(efflist);


features = [IBeff_flat IBlapse_flat IBslope_flat IBoffset_flat;
    Qeff_flat Qlapse_flat Qslope_flat Qoffset_flat];

% Normalize features
features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

[idx, C] = kmeans(features_norm, N);

% Permute the indices for a more logical order
% We want: 2->5, 5->2

idxcopy = idx;
idx(idxcopy == 2) = 5;
idx(idxcopy == 5) = 2;

% Bring back to orginial shape
idxIB = idx(1:numel(IBeff_flat));
idxQ = idx(numel(IBeff_flat) + 1:end);
idxIB =reshape(idxIB, IBxdim, IBydim);
idxQ = reshape(idxQ, Qxdim, Qydim);

% figure;
% subplot(121)
produce_heatmap(idxIB, prewlst, pswitchlst, 'clim', [1 N], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);

% subplot(122)
produce_heatmap(idxQ, epslst, gammalst, 'clim', [1 N], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');



%% SVM on the class labels
% Load simult sim dataset
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/svmresults_from_pickle_083121.mat');
% idxQresize = ceil(imresize(idxQ, [11 11]));
% idxIBresize = floor(imresize(idxIB, [11 11]));

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


features = [IBeffall IBlapseall IBslopeall IBoffsetall;
    Qeffall Qlapseall Qslopeall Qoffsetall];
features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

%shuffle
order = randperm(numel(labels));
labels_shuffled = labels(order);
features_shuffled = features(order,:);


%80% training, 20% testing
rng('shuffle');
ntrain = floor(numel(labels) * 0.8);
Xtrain = features_shuffled(1:ntrain,:);
ytrain = labels_shuffled(1:ntrain);
Xtest = features_shuffled(ntrain + 1:end,:);
ytest = labels_shuffled(ntrain + 1:end);


% t = templateLinear();
t = templateSVM('Standardize',true, 'KernelFunction', 'rbf');
Mdl = fitcecoc(Xtrain',ytrain,'Learners',t,'ObservationsIn','columns');
ypred = Mdl.predict(Xtest);

% Performance
perf = sum(ypred == ytest) / numel(ytest);

% Make the confusion matrix
counts = nan(5,5);
for i = 1:5
    for j = 1:5
        counts(i,j) = sum(ypred == i & ytest == j);

    end
end

confusion = counts ./ sum(counts, 1);


%%
confusionchart(ytest, ypred,'RowSummary','row-normalized');
% mymakeaxis('x_label', 'Predicted Class');
set(gca,'FontSize', 16, 'FontName', 'helvetica')


%% k-means visualizing
% figure;
% for i = 1:4
%     for j = i:4
%         visualize_kmeans(i,j, idx, features_norm)
%         title(sprintf('i = %d, j = %d', i, j));
% %         visualize_kmeans(1,3, idx, features_norm)
% %         visualize_kmeans(1,4, idx, features_norm)
%     end
% end

%% Regime characterization
features1 = features(idx == 1,:);
features2 = features(idx == 2,:);
features3 = features(idx == 3,:);
features4 = features(idx == 4,:);
features5 = features(idx == 5,:);


means = [mean(features1, 1) mean(features2, 1) mean(features3, 1), ...
    mean(features4, 1), mean(features5, 1)];
stds = [std(features1, [], 1), std(features2, [], 1), std(features3, [], 1), ...
    std(features4, [], 1), std(features5, [], 1)];

% figure;
% errorbar(1:5, means(1:4:end), stds(1:4:end), 'o', 'MarkerFaceColor', 'b',...
%     'LineWidth', 2, 'MarkerEdgeColor', 'b', 'Color', 'b');
% mymakeaxis('x_label', 'Class', 'y_label', 'Efficiency')
make_feature_plot(1, means, stds, 'Efficiency', nan);
make_feature_plot(2, means, stds, 'Lapse', [0, 0.3]);
make_feature_plot(3, means, stds, 'Slope', [0, 20]);
make_feature_plot(4, -means, stds, 'Offset', nan);


% figure;
% errorbar



%%

% Decoding of high-gamma regime vs good inf-based regime
gammaid = 10;
epsid = 5;

psid = 8;
prid = 7;

Qeff = squeeze(Qeff_arr(gammaid, epsid, :));
Qlapse = squeeze(Qlapse_arr(gammaid, epsid, :));
Qslope = squeeze(Qslope_arr(gammaid, epsid, :));
Qoffset = squeeze(Qoffset_arr(gammaid, epsid, :));

IBeff = squeeze(IBeff_arr(prid, psid, :));
IBlapse = squeeze(IBlapse_arr(prid, psid, :));
IBslope = squeeze(IBslope_arr(prid, psid, :));
IBoffset = squeeze(IBoffset_arr(prid, psid, :));


figure;
subplot(221)
plot(Qoffset)
hold on
plot(IBoffset);

subplot(222)
plot(Qslope)
hold on
plot(IBslope);

subplot(223)
plot(Qlapse)
hold on
plot(IBlapse);

subplot(224)
plot(Qeff)
hold on
plot(IBeff);






function visualize_kmeans(i, j, idx, features_norm)
subplot(4,4,(i-1) * 4 + j)
plot(features_norm(idx == 1, i), features_norm(idx == 1, j), '.');
hold on
plot(features_norm(idx == 2, i), features_norm(idx == 2, j), '.');
plot(features_norm(idx == 3, i), features_norm(idx == 3, j), '.');

end

function make_feature_plot(classid, means, stds, featurename, ylims)
figure('Position', [440,376,320,422]);
paperaesthetics;
colors = brewermap(10, 'Blues');
% coltouse = colors(8, :);
coltouse = [8, 81, 156]/ 255;
errorbar(1:5, means(classid:4:end), stds(classid:4:end), 'o', 'MarkerFaceColor', coltouse,...
    'LineWidth', 2, 'MarkerEdgeColor', coltouse, 'Color', coltouse);

if ~isnan(ylims)
    ylim(ylims)
end

mymakeaxis('x_label', 'Class', 'y_label', featurename, 'xticks', 1:5)
end



