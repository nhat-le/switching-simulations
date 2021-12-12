%% Script for standardizing classification and segmentation for all probabilities

% Load raw data
opts = struct;
opts.version = '121021';
[res1, opts] = load_and_run(0, opts);
% res2 = load_and_run(0.1, opts);
% res3 = load_and_run(0.2, opts);
% res4 = load_and_run(0.3, opts);


%% Decoding analysis (all probabilities)
opts.reps = 20;
opts.method = 'knn';
opts.nNeighbors = 1;
opts.svmdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/svm/';
opts.svm_version = '121021';
opts.save_model = 1;


[counts_allprob1, Mdls1, MCCs] = do_decoding(1, res1, opts);

% opts.method = 'svm';

%%
MCCs_means = [];
MCCs_stds = [];
Models = {};

for i = 1:11
    if i == 11
        opts.method = 'svm';
    else
        opts.method = 'knn';
        opts.nNeighbors = i;
    end
    [~, Mdl, MCCs] = do_decoding(1, res1, opts);
    MCCs_means(i) = mean(MCCs);
    MCCs_stds(i) = std(MCCs);
    Models{i} = Mdl;
    
end





% [counts_allprob09, Mdls09] = do_decoding(0.9, res2new, opts);
% [counts_allprob08, Mdls08] = do_decoding(0.8, res3new, opts);
% [counts_allprob07, Mdls07] = do_decoding(0.7, res4new, opts);

%%
savedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/svm/models';
filename = sprintf('decoding_common_%s_with%sMdl.mat', opts.svm_version, opts.method);
savename = fullfile(savedir, filename);
if opts.save_model
    if ~exist(savename, 'file')
%         save(savename, 'counts_allprob1', 'counts_allprob09', 'counts_allprob08',...
%             'counts_allprob07', 'Mdls1', 'Mdls09', 'Mdls08', 'Mdls07')
        save(savename, 'counts_allprob1', 'Mdls1')

    else
        error('File exists')
    end
end

function [counts_all, Mdls, MCCs] = do_decoding(prob, res, opts)

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

% Qoffsetall(Qoffsetall < -20) = 3;
% IBoffsetall(IBoffsetall < -20) = 3;



features = [IBeffall IBlapseall IBslopeall IBoffsetall;
    Qeffall Qlapseall Qslopeall Qoffsetall];
% features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

% To balance the number of examples for each class
counts = [];
for ilabel = 1:5
    counts(ilabel) = sum(labels == ilabel);
end

mincounts = min(counts);
filteredIDs = [];
for ilabel = 1:5
    label_pos = find(labels == ilabel);
    filtered_lbl_pos = randsample(label_pos, mincounts, false);
    filteredIDs = [filteredIDs; filtered_lbl_pos];
end


%shuffle
counts_all = {};
Mdls = {};
MCCs = [];

for k = 1:opts.reps
%     order = randperm(numel(labels));
    fprintf('Repetition %d\n', k);
    order = randsample(filteredIDs, numel(labels), true);
    labels_shuffled = labels(order);
    features_shuffled = features(order,:);

    %80% training, 20% testing
    rng('shuffle');
    ntrain = floor(numel(filteredIDs) * 0.8);
    Xtrain = features_shuffled(1:ntrain,:);
    ytrain = labels_shuffled(1:ntrain);
    Xtest = features_shuffled(ntrain + 1:end,:);
    ytest = labels_shuffled(ntrain + 1:end);

    t = templateSVM('Standardize',true, 'KernelFunction', 'rbf');
    
    Xtrainraw = res.features;
    Xtrainraw(:,3) = -Xtrainraw(:,3);
    ytrainraw = res.idx;
    
    % TODO: check if need to transpose for SVM model, if not, do an if else
    % condition check
    
    if strcmp(opts.method, 'knn')
        Mdl = fitcknn(Xtrain,ytrain,'NumNeighbors', opts.nNeighbors,'Standardize',1);
        ypred = Mdl.predict(Xtest);
    else
        Mdl = fitcecoc(Xtrain',ytrain,'Learners',t,'ObservationsIn','columns');
        ypred = Mdl.predict(Xtest);
    end

    % Make the confusion matrix
    N = numel(unique(ypred));
    fprintf('Number of clusters = %d\n', N);
    
    counts = confusionmat(ytest, ypred); %confusion matrix
    MCCs(k) = matthews_corr(counts');
    % Matthews correlation coefficient
    
    
    counts_all{k} = counts';
    Mdls{k} = Mdl;
    
end
end