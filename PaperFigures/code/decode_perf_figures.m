% Figures for decoding performance
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/svmresults_from_pickle_090221.mat');
% 
% mindecodingQ = min(decoding_perf, [], [3, 4]);
% mindecodingIB = squeeze(min(decoding_perf, [], [1, 2]));
% 
% %%
% produce_heatmap(mindecodingQ', epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Decoding performance', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');
% 
% 
% %%
% produce_heatmap(mindecodingIB', prewlst, pswitchlst, 'clim', [0.5, 1], 'legendname', 'Decoding performance', ...
%     'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);


%% Regime demarcation (with k-means clustering)
% Form the feature vectors
% filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat';
addpath('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code')
prob = 0.7;
rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata';
% filedir = sprintf('%s/EGreedyQLearningAgent-withCorr-prob%.2fto%.2f-072321.mat', rootdir, 1-prob, prob);
filedir = sprintf('%s/EGreedyQLearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-092321.mat', rootdir, 1-prob, prob);

load(filedir);

% rng('shuffle')
rng(127)
N = 5;
Qeff_flat = reshape(efflist, [], 1);
Qlapse_flat = reshape(LapseL, [], 1);
Qslope_flat = reshape(PLslopelist, [], 1);
Qoffset_flat = reshape(PLoffsetlist, [], 1);


% Remove outliers from slope (commonly seen)
% upperQslope = median(Qslope_flat) + 2 * std(Qslope_flat);
% Qslope_flat(Qslope_flat > upperQslope) = nan;

[Qxdim, Qydim] = size(efflist);
% filedir = sprintf('%s/EGreedyInferenceBasedAgent-withCorr-prob%.2fto%.2f-072121.mat', rootdir, 1-prob, prob);
filedir = sprintf('%s/EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-092321.mat', rootdir, 1-prob, prob);

% filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-withCorr-prob0.00to1.00-072121.mat';
load(filedir);

IBeff_flat = reshape(efflist, [], 1);
IBlapse_flat = reshape(LapseL, [], 1);
IBslope_flat = reshape(PLslopelist, [], 1);
IBoffset_flat = reshape(PLoffsetlist, [], 1);

IBoffset_flat(IBoffset_flat < -20) = nan;
Qoffset_flat(Qoffset_flat < -20) = nan;

% for 0.7 prob
% upperIBslope = median(IBslope_flat) + 1.6 * std(IBslope_flat);
upperIBslope = median(IBslope_flat) + 2 * std(IBslope_flat);

IBslope_flat(IBslope_flat > upperIBslope) = nan;


[IBxdim, IBydim] = size(efflist);


features = [IBeff_flat IBlapse_flat IBslope_flat IBoffset_flat;
    Qeff_flat Qlapse_flat Qslope_flat Qoffset_flat];

% Normalize features
features_norm = (features - nanmean(features, 1)) ./ nanstd(features, [], 1);

[idx, C] = kmeans(features_norm, N, 'MaxIter', 2000);


% Refine the labels..
idx(idx == 3 & features(:,4) > -2) = 4;


idxcopy = idx;

%for re-arranging labels
switch prob
    case 1
%         idx(idxcopy == 2) = 5;
%         idx(idxcopy == 5) = 2;
%         idx = rotate(idx, [2, 3]);
    case 0.9
%         idx(idxcopy == 3) = 2;
%         idx(idxcopy == 5) = 3;
%         idx(idxcopy == 2) = 5;
%         idx = rotate(idx, [3, 2, 5]);
    case 0.8
%         idx = rotate(idx, [3, 1, 2, 5]);
    case 0.7
%         idx = rotate(idx, [1, 5]);
%         idx = rotate(idx, [3, 4]);
end

% Bring back to orginial shape
idxIB = idx(1:numel(IBeff_flat));
idxQ = idx(numel(IBeff_flat) + 1:end);
idxIB =reshape(idxIB, IBxdim, IBydim);
idxQ = reshape(idxQ, Qxdim, Qydim);





% idxIB(isnan(idxIB)) =  N+1;
% idxQ(isnan(idxQ)) = N+1;


% figure;
% subplot(121)
cmap = brewermap(N, 'Blues');
cmap = [cmap; 0.2 0.2 0.2];
h = produce_heatmap(idxIB, prewlst, pswitchlst, 'clim', [0.5 N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, N + 0.5]);

% Save Q plot
filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs/perfRegimeIB_prob%.1f.pdf', 1-prob);
if ~exist(filename, 'file')
    saveas(gcf, filename);
end


% subplot(122)
produce_heatmap(idxQ, epslst, gammalst, 'clim', [0.5 N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, N + 0.5]);
% Save IB plot
filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs/perfRegimeQ_prob%.1f.pdf', 1-prob);
if ~exist(filename, 'file')
    saveas(gcf, filename);
end


% Visualize clusters
figure;
hold on
for i = 1:4
    for j = 1:4
        subplot(4,4,(i-1)*4+j)
        hold on
        plot(features(idx==1,i), features(idx==1,j), '.'); %blue
        plot(features(idx==2,i), features(idx==2,j), '.'); %red
        plot(features(idx==3,i), features(idx==3,j), '.'); %yellow
        plot(features(idx==4,i), features(idx==4,j), '.'); %purple
        plot(features(idx==5,i), features(idx==5,j), '.'); %green
%         plot(features(idx==1 & features(:,4) > -2,i), features(idx==1 & features(:,4) > -2,j), 'kx'); %green
        
        if i == 1
            xlabel('Eff')
        elseif i == 2
            xlabel('Lapse')
        elseif i == 3
            xlabel('Slope')
        elseif i == 4 
            xlabel('Offset')
        end
        
    end
end


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


%% SVM on the class labels
% Load simult sim dataset
load(sprintf('%s/svmresults_from_pickle_092221_prob%.2f.mat', rootdir, 1-prob));
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

Qoffsetall(Qoffsetall < -20) = -20;
IBoffsetall(IBoffsetall < -20) = -20;



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

% save('svm_Mdl.mat', 'Mdl', 'idxQ' , 'idxIB');

%% Visualize clusters and labels
figure;
hold on
for i = 1:4
    for j = 1:4
        subplot(4,4,(i-1)*4+j)
        hold on
        plot(features(labels==1,i), features(labels==1,j), '.'); %blue
        plot(features(labels==2,i), features(labels==2,j), '.'); %red
        plot(features(labels==3,i), features(labels==3,j), '.'); %yellow
        plot(features(labels==4,i), features(labels==4,j), '.'); %purple
        plot(features(labels==5,i), features(labels==5,j), '.'); %green

        if i == 1
            xlabel('Eff')
        elseif i == 2
            xlabel('Lapse')
        elseif i == 3
            xlabel('Slope')
        elseif i == 4 
            xlabel('Offset')
        end
        
    end
end


%%
figure;
hold off
confusionchart(ytest, ypred,'RowSummary','row-normalized');
% mymakeaxis('x_label', 'Predicted Class');
set(gca,'FontSize', 16, 'FontName', 'helvetica')









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


function res = rotate(arr, order)
res = arr;
for i = 1:numel(order) - 1
    res(arr == order(i)) = order(i+1);
end

res(arr == order(end)) = order(1);

end

