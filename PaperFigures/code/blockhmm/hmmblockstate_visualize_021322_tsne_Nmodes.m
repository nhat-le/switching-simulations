% script for official paper figures sent on 2.21.2022

paths = pathsetup('matchingsim');
expfitdate = '121821bK3';
opts.rootdir = fullfile(paths.blockhmmfitpath, expfitdate);
folders = dir(fullfile(opts.rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
mdltypes = 1:2;
mdlids = 1:10;
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.python_assist = 0;
opts.effmethod = 'sim';
opts.savefile = 0;
opts.showplot = 0;
% opts.model_name = 'decoding_common_010522_withsvmMdl_knn_svm_v10_tsne.mat';
opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10_tsne.mat';
opts.svmmodelpath = paths.svmmodelpath;


[~, aggmeans_native, aggparams_native] = load_params(folders, opts);

% concat and do k-means
aggparams_all = cell2mat(aggparams_native)';
aggmeans_all = cell2mat(aggmeans_native');

yvals_all = [];
for k = 1:size(aggparams_all, 1)
    param = aggparams_all(k,:);
    yvals_all(k,:) = mathfuncs.sigmoid(1:15, param(1), param(2), param(3));
%         plot(yvals, 'b');
end


%% clustering and plotting of the HMM modes
seed = 19;
rng(seed);
K = 6;
idx = kmeans(aggmeans_all, K);
%shuffle the idx
idxcopy = idx;
% idx(idxcopy == 6) = 2;
% idx(idxcopy == 2) = 6;
% 
% idx(idxcopy == 4) = 3;
% idx(idxcopy == 5) = 4;
% idx(idxcopy == 3) = 5;

% idx(idxcopy == 3) = 6;
% idx(idxcopy == 4) = 2;
% idx(idxcopy == 5) = 1;
% idx(idxcopy == 6) = 4;



figure;
% plotting raw
for i = 1:K   
    aggmeans_sub = yvals_all(idx == i,:);
    meanall = mean(aggmeans_sub(:));
    
%     figure;
    subplot(2,3,i);
    plot(aggmeans_sub', 'Color', 'k', 'LineWidth', 0.25);
    hold on
    plot(mean(aggmeans_sub, 1), 'r', 'Color', 'r', 'LineWidth', 2);
%     mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)');
%     saveas(gcf, fullfile(paths.figpath, sprintf('hmmblockFigs/hmmblock_121821_mode%d.pdf', i)));

end


%% plotting smooth
% plotting raw
for i = 1:K   
    param_sub = aggparams_all(idx == i,:);
    yvals_sub = [];
    for k = 1:size(param_sub, 1)
        param = param_sub(k,:);
        yvals_sub(k,:) = mathfuncs.sigmoid(1:15, param(1), param(2), param(3));
%         plot(yvals, 'b');
    end
    
    
%     figure;
    subplot(2,3,i);
    plot(yvals_sub', 'k', 'LineWidth', 0.25)
    hold on
    plot(mean(yvals_sub, 1), 'r', 'Color', 'r', 'LineWidth', 2);
    ylim([0, 1])
    mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)', ...
        'yticks', 0:0.2:1, 'xticks', 0:5:15);


end



%% Decoding into inf-based/qlearning strategies
load(fullfile(opts.svmmodelpath, opts.model_name), 'Models');

num_state5 = [];
for i = 24
    opts.mdltype = i;
    opts.mdlid = 5;
    [all_aggparams, aggmeans_all, statesFlat, features_flat] = apply_model(aggparams_native, aggmeans_native, opts);

    % Plot
    fprintf('i=%d, num states 6 = %d\n', i, sum(statesFlat == 6));
end

% Plot state composition
composition = nan(max(idx), max(statesFlat));
for i = 1:max(idx)
    for j = 1:max(statesFlat)
        composition(i,j) = sum(idx == i & statesFlat == j);  
    end 
end


%% Visualize the transition functions grouped by strategy
figure;
for i = 1:6
    subplot(2,3,i)
    mode_funcs = yvals_all(statesFlat == i,:)';
    plot(mode_funcs, 'k', 'LineWidth', 0.25);
    hold on
    plot(mean(mode_funcs, 2), 'r', 'Color', 'r', 'LineWidth', 2);

    ylim([0, 1])
    mymakeaxis('x_label', 'Trials in block', 'y_label', 'P(Correct)', ...
        'yticks', 0:0.2:1, 'xticks', 0:5:15);
end



%%
statefrac = composition ./ sum(composition, 2);
cols = paperaesthetics;

if opts.showplot
    figure('Position', [440,423,736,375]);
    h = bar(statefrac,'stacked');
    for i = 1:6
        h(i).FaceColor = cols.colors(i,:);
        h(i).ShowBaseLine = 'off';
    end
    % xlim([0.5 4.5])
    mymakeaxis('x_label', 'HMM mode', 'y_label', 'Fraction', 'xticks', 1:6)
    l = legend(h(1:6), {'Q-learning 1', 'Q-learning 2', 'Q-learning 3', 'Q-learning 4', 'Inference-based 5',...
        'Inference-based 6'}, 'Position', [0.4,0.42,0.1,0.1], ...
        'FontSize', 14);
end


%% Save if requested
savefilename = fullfile(opts.rootdir, sprintf('hmmblock_classification_info_%s_v1.mat', expfitdate));
if opts.savefile && ~exist(savefilename, 'file')
    blockhmm_idx = idx;
    save(savefilename, 'opts', 'statesFlat', 'blockhmm_idx', 'folders', 'aggmeans_native', 'aggparams_native');
    fprintf('file saved!\n')
else
    fprintf('file exists, skipping save...\n');
end

    
    