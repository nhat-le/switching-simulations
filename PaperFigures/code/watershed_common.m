%% Script for standardizing classification and segmentation for all probabilities

% Load raw data
[res1, opts] = load_and_run(0);
res2 = load_and_run(0.1);
res3 = load_and_run(0.2);
res4 = load_and_run(0.3);

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


% Characterize each cluster for each prob
meansAll = {};
stdsAll = {};
Nclust = numel(unique(res1new.idx));
Nprobs = 4;
res_all = {res1new, res2new, res3new, res4new};
for i = 1:Nprobs
    means = [];
    stds = [];
    for nclust = 1:Nclust
        featuresi = res_all{i}.features(res_all{i}.idx == nclust,:);
        means(nclust,:) = mean(featuresi, 1);
        stds(nclust,:) = std(featuresi, [], 1);
        
    end
    
    meansAll{i} = means;
    stdsAll{i} = stds;
end


%% MDS plot

for probi = 1:4
    figure;
    hold on
    for i = 1:Nclust
        col2use = colors(i,:);
        if probi == 2
            plot(res_all{probi}.Y(res_all{probi}.idx == i, 1), -res_all{probi}.Y(res_all{probi}.idx == i, 2), 'o', ...
                'MarkerFaceColor', col2use, 'MarkerEdgeColor', 'w', 'MarkerSize', 7)            
        else
            plot(res_all{probi}.Y(res_all{probi}.idx == i, 1), res_all{probi}.Y(res_all{probi}.idx == i, 2), 'o', ...
                'MarkerFaceColor', col2use, 'MarkerEdgeColor', 'w', 'MarkerSize', 7)
        end
    end
    mymakeaxis('x_label', 'PC1', 'y_label', 'PC2')
end

%% IB/Q space heat plots






%% Plotting
figure('Position', [626,401,721,353]);
hold on
paperaesthetics;
% colors = brewermap(6, 'BuGn');
colors = brewermap(6, 'Set1');
xvals = linspace(-0.15, 0.15, Nclust);

plottype = {'Efficiency', 'Lapse', 'Slope', 'Offset'};
k = 4;

handles = [];
for i = 1:numel(meansAll)
    coltouse = colors(i,:);
    
    means = meansAll{i};
    if k == 4
        means = -means;
    end
    stds = stdsAll{i};
    h = errorbar((1:Nclust) + xvals(i)', means(:,k), stds(:,k), 'o', 'MarkerFaceColor', coltouse,...
    'LineWidth', 0.75, 'MarkerEdgeColor', coltouse, 'Color', coltouse);
    handles(i) = h;
end

switch k
    case 1
        ylim([0.5 1])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:Nclust, 'yticks', 0.5:0.1:1)

    case 2
        ylim([0 0.5])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:Nclust, 'yticks', 0:0.1:0.5)

    case 3
        ylim([0 3.3])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:Nclust, 'yticks', 0:1:3)
    
    case 4
        ylim([0 14])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:Nclust)
    
        
end
        
if k < 4
    l = legend(handles, {'1', '0.9', '0.8', '0.7'}, 'Position', [0.3766 0.7337 0.0749 0.1827]);
    l.Title.String = 'Probability';
    l.Color = 'none';
else
    l = legend(handles, {'1', '0.9', '0.8', '0.7'});
    l.Title.String = 'Probability';
    l.Color = 'none';
    
end





%% Decoding, can skip after processing
% Load simult sim dataset
opts.prob = 1;
load(sprintf('%s/svmresults_from_pickle_092221_prob%.2f.mat', opts.rootdir, 1-opts.prob));

[idxQ, idxIB] = reshapeidx(idxprob1, res1new);


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



features = [IBeffall IBlapseall IBslopeall IBoffsetall;
    Qeffall Qlapseall Qslopeall Qoffsetall];
features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

%shuffle
confusion_all = {};
for k = 1:20
%     order = randperm(numel(labels));
    fprintf('Repetition %d\n', k);
    order = randsample(1:numel(labels), numel(labels), true);
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
    N = numel(unique(ypred));
    counts = nan(N,N);
    for i = 1:N
        for j = 1:N
            counts(i,j) = sum(ypred == i & ytest == j);

        end
    end

    counts_all{k} = counts;
    
end

% save
% save('prob1_counts.mat', 'counts_all');


%% Now post-process
counts_allprob1 = do_decoding(1, res1new, opts);
counts_allprob09 = do_decoding(0.9, res2new, opts);
counts_allprob08 = do_decoding(0.8, res3new, opts);
counts_allprob07 = do_decoding(0.7, res4new, opts);

%% Plot
load('decoding_common_092821.mat');
counts_all = {counts_allprob1, counts_allprob09, counts_allprob08, counts_allprob07};
% colors = brewermap(6, 'BuGn');
colors = brewermap(6, 'Set1');
figure;
hold on
for probi = 1:numel(counts_all)
    allprob = counts_all{probi};
    means = [];
    stds = [];
    coltouse = colors(probi,:);
    Nclust = size(allprob{1}, 1);
    for i = 1:Nclust
        perf_clusti = cellfun(@(x) find_perf(x, i), allprob);
        means(i) = mean(perf_clusti);
        stds(i) = std(perf_clusti);
    end


    h = errorbar(1:Nclust, means, stds, 'o', 'LineWidth', 0.75, ...
        'MarkerEdgeColor', coltouse, 'Color', coltouse);
    handles(probi) = h;
    plot(1:Nclust, means, 'Color', coltouse, 'LineWidth', 0.75);
%     plot(mapvals{probi}, means, 'LineWidth', 2);
end

mymakeaxis('x_label', 'Class', 'y_label', 'Decoding accuracy')
l = legend(handles, {'1', '0.9', '0.8', '0.7'}, 'Position', [0.46,0.41,0.12,0.19]);
l.FontSize = 12;
l.Title.String = 'Probability';
l.Title.FontSize = 12;
l.Color = 'none';



function counts_all = do_decoding(prob, res, opts)
opts.prob = prob;
load(sprintf('%s/svmresults_from_pickle_092221_prob%.2f.mat', opts.rootdir, 1-opts.prob));

[idxQ, idxIB] = reshapeidx(res.idx, res);


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



features = [IBeffall IBlapseall IBslopeall IBoffsetall;
    Qeffall Qlapseall Qslopeall Qoffsetall];
% features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

%shuffle
counts_all = {};
for k = 1:20
%     order = randperm(numel(labels));
    fprintf('Repetition %d\n', k);
    order = randsample(1:numel(labels), numel(labels), true);
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
    N = numel(unique(ypred));
    fprintf('Number of clusters = %d\n', N);
    counts = nan(N,N);
    for i = 1:N
        for j = 1:N
            counts(i,j) = sum(ypred == i & ytest == j);

        end
    end

    counts_all{k} = counts;
    
end
end

function res = find_perf(arr, i)
res = arr(i,i) / sum(arr(:,i));

end


function [out, opts] = load_and_run(prob)
folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';

switch prob
    case 0
        filename = 'opts_prob0.0-2021-09-25 20.52.mat';
    case 0.1
        filename = 'opts_prob0.1-2021-09-25 21.44.mat';
    case 0.2
        filename = 'opts_prob0.2-2021-09-25 21.57.mat';
    case 0.3
        filename = 'opts_prob0.3-2021-09-25 22.29.mat';
end

load(fullfile(folder, filename));
opts.save = 0;
opts.savefeatures = 0;


[idx, out] = run_watershed(opts);

% Rotation for idx
switch prob
    case 0
        idx = rotate(idx, [3, 2, 4]);
    case 0.1
        idx = rotate(idx, [6, 1, 3]);
        idx = rotate(idx, [4 5]);
    case 0.2
        idx = rotate(idx, [2, 1, 4]);
        idx = rotate(idx, [5, 3]);
    case 0.3
        idx = rotate(idx, [4, 1]);
        idx = rotate(idx, [2, 5, 3]);
end

%[idxQ, idxIB] = reshapeidx(idx, out);
out.idx = idx;


end





function res = rotate(arr, order)
res = arr;
for i = 1:numel(order) - 1
    res(arr == order(i)) = order(i+1);
end

res(arr == order(end)) = order(1);

end