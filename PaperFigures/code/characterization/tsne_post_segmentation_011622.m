%% Script for standardizing classification and segmentation for all probabilities
% Load raw data
paths = pathsetup('matchingsim');
% Load raw data
opts = struct;
opts.version = '092321';

[res1, opts_out] = load_and_run_tsne(0, opts);
res2 = load_and_run_tsne(0.1, opts);
res3 = load_and_run_tsne(0.2, opts);
res4 = load_and_run_tsne(0.3, opts);



%%
% Characterize each cluster for each prob
meansAll = {};
stdsAll = {};
Nclust = max(res1.idx);
res_all = {res1, res2, res3, res4};
for i = 1:numel(res_all)
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


%% MDS plot without colors

% figure()
%     
% hold on
% for i = 1:Nclust
% %     col2use = colors(i,:);
%     plot(res1.Y(:, 1), res1.Y(:, 2), 'o', ...
%         'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'MarkerSize', 7)             
% end
% mymakeaxis('x_label', 'PC1', 'y_label', 'PC2')

%% MDS plot colored by cluster
colors = brewermap(6, 'Set1');
probs = [1, 0.9, 0.8, 0.7];
colors = colors([2,1,5,6,4,3],:); %permute to align with MATLAB default..
for probi = 1:4
    figure()
    
    hold on
    for i = 1:Nclust
        col2use = colors(i,:);
        if probi == 2 || probi == 1
            plot(res_all{probi}.Y(res_all{probi}.idx == i, 1), -res_all{probi}.Y(res_all{probi}.idx == i, 2), 'o', ...
                'MarkerFaceColor', col2use, 'MarkerEdgeColor', 'w', 'MarkerSize', 7)            
        else
            plot(res_all{probi}.Y(res_all{probi}.idx == i, 1), res_all{probi}.Y(res_all{probi}.idx == i, 2), 'o', ...
                'MarkerFaceColor', col2use, 'MarkerEdgeColor', 'w', 'MarkerSize', 7)
        end
    end
    mymakeaxis('x_label', 'PC1', 'y_label', 'PC2')
end

%% IB/Q space regime demarcation
savefig = 0;
currdate = datetime;
currdate.Format = 'yyyy-MM-dd HH.mm';
currdate = string(currdate); 
for i = 1:4
    resstruct = res_all{i};
    [idxQ, idxIB] = reshapeidx(resstruct.idx, resstruct);
    prewlst = resstruct.prewlst;
    pswitchlst = resstruct.pswitchlst;
    epslst = resstruct.epslst;
    gammalst = resstruct.gammalst;
    % figure;
    % subplot(121)
    cmap = brewermap(6, 'Set1');
%     cmap = brewermap(12, 'Paired');
%     cmap = cmap(1:2:end,:);
%     cmap = cmap([1, 3, 4, 5, 2, 6],:); %permute to align with MATLAB default..
    % cmap = [cmap; 0.2 0.2 0.2];
    h = produce_heatmap(idxIB, prewlst, pswitchlst, 'clim', [0.5 Nclust+1.5], 'legendname', 'Performance regime', ...
        'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, Nclust + 0.5]);
    filename = sprintf('%s/perfRegimeIB_common_prob%.1f-%s.pdf', opts.savepath, 1-probs(i), currdate);
%     if ~exist(filename, 'file')
%         saveas(gcf, filename);
%     end

    produce_heatmap(idxQ, epslst, gammalst, 'clim', [0.5 Nclust+1.5], 'legendname', 'Performance regime', ...
        'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, Nclust + 0.5]);
    filename = sprintf('%s/perfRegimeQ_common_prob%.1f-%s.pdf', opts.savepath, 1-probs(i), currdate);
%     if ~exist(filename, 'file')
%         saveas(gcf, filename);
%     end

end



%% Feature characterization
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
    h = errorbar((1:Nclust) + xvals(i)', means(:,k), stds(:,k), 'o-', 'MarkerFaceColor', coltouse,...
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

%% Feature characterization v2 (grouped by clusters)
probs = [0 0.1 0.2 0.3];
figure('Position', [626,401,721,353]);
hold on
paperaesthetics;
cols = paperaesthetics;
colors = cols.colors; 

Nprobs = numel(meansAll);

xvals = linspace(-0.015, 0.015, Nclust);

plottype = {'Efficiency', 'Lapse', 'Slope', 'Offset'};
k = 4;

handles = [];

means = [];
stds = [];
for j = 1:Nprobs
    means = [means meansAll{j}(:,k)];
    stds = [stds stdsAll{j}(:,k)];
end
if k == 4
    means = -means;
end

for i = 1:Nclust
    coltouse = colors(i,:);
  
    if i == 4
        h = errorbar(probs + xvals(i)', means(i,:), stds(i,:), 'o-', 'MarkerFaceColor', coltouse,...
        'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'Color', 'k');
    else
        h = errorbar(probs + xvals(i)', means(i,:), stds(i,:), 'o-', 'MarkerFaceColor', coltouse,...
        'LineWidth', 0.75, 'MarkerEdgeColor', coltouse, 'Color', coltouse);
    end
    handles(i) = h;
end

switch k
    case 1
        plot(probs, [1, 0.9, 0.8, 0.7], 'k--')
        plot(probs, [0.5 0.5 0.5 0.5], 'k--')
        
        ylim([0.5 1])
        mymakeaxis('x_label', 'Probability', 'y_label', plottype{k}, 'xticks', probs, 'yticks', 0.5:0.1:1,...
            'xticklabels', {'100-0', '90-10', '80-20', '70-30'}, 'font_size', 22)

    case 2
        ylim([0 0.5])
        mymakeaxis('x_label', 'Probability', 'y_label', plottype{k}, 'xticks', probs, 'yticks', 0:0.1:0.5,...
            'xticklabels', {'100-0', '90-10', '80-20', '70-30'}, 'font_size', 22)

    case 3
        ylim([0 3.3])
        mymakeaxis('x_label', 'Probability', 'y_label', plottype{k}, 'xticks', probs, 'yticks', 0:1:3,...
            'xticklabels', {'100-0', '90-10', '80-20', '70-30'}, 'font_size', 22)
    
    case 4
        ylim([0 14])
        mymakeaxis('x_label', 'Probability', 'y_label', plottype{k}, 'xticks', probs,...
            'xticklabels', {'100-0', '90-10', '80-20', '70-30'}, 'font_size', 22)
        
end


%% Decoding analysis (all probabilities)
opts.reps = 20;
opts.method = 'knn';
opts.nNeighbors = 5;
[counts_allprob1, Mdls1] = do_decoding(1, res1new, opts);
[counts_allprob09, Mdls09] = do_decoding(0.9, res2new, opts);
[counts_allprob08, Mdls08] = do_decoding(0.8, res3new, opts);
[counts_allprob07, Mdls07] = do_decoding(0.7, res4new, opts);

savedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata';
filename = 'decoding_common_101421_withknnMdl.mat';
savename = fullfile(savedir, filename);
% if ~exist(savename, 'file')
%     save(savename, 'counts_allprob1', 'counts_allprob09', 'counts_allprob08',...
%         'counts_allprob07', 'Mdls1', 'Mdls09', 'Mdls08', 'Mdls07')
% end
%% Plot decoding accuracy, grouped by class
% load('decoding_common_092821.mat');
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
        perf_clusti = cellfun(@(x) find_sensitivity(x, i), allprob);
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

%% Plot decoding accuracy, grouped by prob
% load('decoding_common_092821.mat');
counts_all = {counts_allprob1, counts_allprob09, counts_allprob08, counts_allprob07};
% colors = brewermap(6, 'BuGn');
colors = brewermap(6, 'Set1');
colors = colors([2,1,5,4,3],:); %permute to align with MATLAB default..


Nprobs = numel(counts_all);
Nclust = size(counts_all{1}{1}, 1);
figure;
hold on


for k = 1:Nclust
%     allprob = counts_all{probi};
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


    h = errorbar(0.7:0.1:1, means, stds, 'o-', 'LineWidth', 0.75, ...
        'MarkerEdgeColor', coltouse, 'Color', coltouse);
    handles(k) = h;
%     plot(1:Nclust, means, 'Color', coltouse, 'LineWidth', 0.75);
%     plot(mapvals{probi}, means, 'LineWidth', 2);
end
ylim([0.5, 1])

mymakeaxis('x_label', 'Probability', 'y_label', 'Precision', 'xticks', [0.7 0.8 0.9 1.0],...
            'xticklabels', {'100-0', '90-10', '80-20', '70-30'})% l = legend(handles, {'1', '0.9', '0.8', '0.7'}, 'Position', [0.46,0.41,0.12,0.19]);
l = legend(handles, {'1', '2', '3', '4', '5'}, 'Position', [0.46,0.41,0.12,0.19]);
l.FontSize = 12;
l.Title.String = 'Class';
l.Title.FontSize = 12;
l.Color = 'none';




function [counts_all, Mdls] = do_decoding(prob, res, opts)

if ~isfield(opts, 'method'); opts.method = 'knn'; end
if ~isfield(opts, 'nNeighbors'); opts.nNeighbors = 5; end


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

    % t = templateLinear();
    t = templateSVM('Standardize',true, 'KernelFunction', 'rbf');
    
    Xtrainraw = res.features;
    Xtrainraw(:,3) = -Xtrainraw(:,3);
    ytrainraw = res.idx;
    
    % TODO: check if need to transpose for SVM model, if not, do an if else
    % condition check
    
    if strcmp(opts.method, 'knn')
        Mdl = fitcknn(Xtrain,ytrain,'NumNeighbors', opts.nNeighbors,'Standardize',1);
%         Mdl = fitcknn(Xtrainraw, ytrainraw, 'NumNeighbors', opts.nNeighbors,'Standardize',1);
        ypred = Mdl.predict(Xtest);
    else
        Mdl = fitcecoc(Xtrain',ytrain,'Learners',t,'ObservationsIn','columns');
        ypred = Mdl.predict(Xtest);
    end

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
    Mdls{k} = Mdl;
    
end
end

function res = find_sensitivity(arr, i)
res = arr(i,i) / sum(arr(:,i));

end


function res = find_precision(arr, i)
res = arr(i,i) / sum(arr(i,:));

end