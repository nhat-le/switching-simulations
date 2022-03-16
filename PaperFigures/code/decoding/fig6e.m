%% Script for standardizing classification and segmentation for all probabilities

paths = pathsetup('matchingsim');
% Load raw data
opts = struct;
opts.version = '092321';

[res1, opts_out] = load_and_run_tsne(0, opts);
res2 = load_and_run_tsne(0.1, opts);
res3 = load_and_run_tsne(0.2, opts);
res4 = load_and_run_tsne(0.3, opts);


%% Decoding analysis (all probabilities)
paths = pathsetup('matchingsim');
opts.nClasses = max(res1.idx);
opts.reps = 20;
opts.method = 'knn';
opts.nNeighbors = 1;
opts.svmdir = paths.svmdatapath;
opts.svm_version = '010522';
opts.save_model = 1;

opts.method = 'knn';
opts.nNeighbors = 24;
[confusions10, Mdl10, ~] = do_decoding(1, res1, opts);
[confusions09, Mdl09, ~] = do_decoding(0.9, res2, opts);
[confusions08, Mdl08, ~] = do_decoding(0.8, res3, opts);
[confusions07, Mdl07, ~] = do_decoding(0.7, res4, opts);


%% Plot decoding accuracy, grouped by prob
counts_all = {confusions10, confusions09, confusions08, confusions07};
cols = paperaesthetics;
colors = cols.colors;

Nprobs = numel(counts_all);
Nclust = size(counts_all{1}{1}, 1);
figure;
hold on


for k = 1:Nclust
    means = [];
    stds = [];
    coltouse = colors(k,:);
    
    means = [];
    stds = [];
    for i = 1:Nprobs
        perf_clusti = cellfun(@(x) find_sensitivity(x, k), counts_all{i});
        means(i) = mean(perf_clusti);
        stds(i) = std(perf_clusti);
    end

    if k == 4
        h = errorbar(0.7:0.1:1, means, stds, 'o-', 'LineWidth', 0.75, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', coltouse, 'Color', 'k');
    else
        h = errorbar(0.7:0.1:1, means, stds, 'o-', 'LineWidth', 0.75, ...
            'MarkerEdgeColor', coltouse, 'MarkerFaceColor', coltouse, 'Color', coltouse);
    end
    handles(k) = h;
%     plot(1:Nclust, means, 'Color', coltouse, 'LineWidth', 0.75);
%     plot(mapvals{probi}, means, 'LineWidth', 2);
end

hline(1/6, 'k--');

ylim([0, 1])

mymakeaxis('x_label', 'Probability', 'y_label', 'Decoding accuracy', 'xticks', [0.7 0.8 0.9 1.0],...
            'xticklabels', {'100-0', '90-10', '80-20', '70-30'}, 'font_size', 20)% l = legend(handles, {'1', '0.9', '0.8', '0.7'}, 'Position', [0.46,0.41,0.12,0.19]);
l = legend(handles, {'Q1', 'Q2', 'Q3', 'Q4', 'IB5', 'IB6'}, 'Position', [0.46,0.41,0.12,0.19]);
l.FontSize = 12;
l.Title.String = 'Regime';
l.Title.FontSize = 12;
l.Color = 'none';


%%
savedir = paths.svmmodelpath;
filename = sprintf('decoding_common_%s_with%sMdl_knn_svm_v10b_tsne.mat', opts.svm_version, opts.method);
notes = 'knn models with k = 2:2:30';
savename = fullfile(savedir, filename);
if opts.save_model
    if ~exist(savename, 'file')
        save(savename, 'Mdl10', 'Mdl09', 'Mdl08', 'Mdl07', 'confusions10', ...
            'confusions09', 'confusions08', 'confusions07', 'opts');
        fprintf('File saved!\n');
    else
        error('File exists')
    end
end

function [confusions, Mdls, MCCs] = do_decoding(prob, res, opts)

if ~isfield(opts, 'method'); opts.method = 'knn'; end
if ~isfield(opts, 'nNeighbors'); opts.nNeighbors = 5; end


opts.prob = prob;

load(sprintf('%s/%s/svmresults_from_pickle_%s_prob%.2f.mat', opts.svmdir,...
    opts.svm_version, opts.svm_version, 1-opts.prob));

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

Qoffsetall(Qoffsetall < -20) = -20;
IBoffsetall(IBoffsetall < -20) = -20;


features = [IBoffsetall IBslopeall IBlapseall IBeffall;
    Qoffsetall Qslopeall Qlapseall Qeffall];
% features = [IBeffall IBlapseall IBslopeall IBoffsetall;
%     Qeffall Qlapseall Qslopeall Qoffsetall];
% features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

% To balance the number of examples for each class
counts = [];
for ilabel = 1:opts.nClasses
    counts(ilabel) = sum(labels == ilabel);
end

mincounts = min(counts);
filteredIDs = [];
for ilabel = 1:opts.nClasses
    label_pos = find(labels == ilabel);
    filtered_lbl_pos = randsample(label_pos, mincounts, false);
    filteredIDs = [filteredIDs; filtered_lbl_pos];
end


%shuffle
confusions = {};
Mdls = {};
MCCs = [];

for k = 1:opts.reps
%     order = randperm(numel(labels));
%     fprintf('Repetition %d\n', k);
    order = randsample(filteredIDs, numel(filteredIDs), false);
    labels_shuffled = labels(order);
    features_shuffled = features(order,:);

    %80% training, 20% testing
%     rng('shuffle');
    ntrain = floor(numel(filteredIDs) * 0.8);
    Xtrain = features_shuffled(1:ntrain,:);
    ytrain = labels_shuffled(1:ntrain);
    Xtest = features_shuffled(ntrain + 1:end,:);
    ytest = labels_shuffled(ntrain + 1:end);
    
    % TODO: check if need to transpose for SVM model, if not, do an if else
    % condition check
    
    if strcmp(opts.method, 'knn')
        Mdl = fitcknn(Xtrain,ytrain,'NumNeighbors', opts.nNeighbors,'Standardize',1);
        ypred = Mdl.predict(Xtest);
    else
%         t = templateSVM('Standardize',true, 'KernelFunction', 'rbf');

        t = templateSVM('Standardize',true,'KernelFunction','rbf', 'Standardize', 1);
        Mdl = fitcecoc(Xtrain',ytrain,'Learners',t,'ObservationsIn','columns');
        ypred = Mdl.predict(Xtest);
    end

    % Make the confusion matrix
    N = numel(unique(ypred));
%     fprintf('Number of clusters = %d\n', N);
    
    counts = confusionmat(ytest, ypred); %confusion matrix
    MCCs(k) = matthews_corr(counts');
    % Matthews correlation coefficient
    
    
    confusions{k} = counts';
    Mdls{k} = Mdl;
    
end
end


function res = find_sensitivity(arr, i)
res = arr(i,i) / sum(arr(:,i));

end


function res = find_precision(arr, i)
res = arr(i,i) / sum(arr(i,:));

end