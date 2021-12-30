path = pathsetup('matchingsim');

prob = 0;

rootdir = fullfile(path.simdatapath, sprintf('122821/EGreedyinf-based*%.2f*.mat', prob));
files = dir(rootdir);
load(fullfile(files(1).folder, files(1).name));
IB_coefs = reshape(coefs_all, size(coefs_all, 1), size(coefs_all, 2), [], 3);

rootdir = fullfile(path.simdatapath, sprintf('122821/EGreedyqlearning*%.2f*.mat', prob));
files = dir(rootdir);
load(fullfile(files(1).folder, files(1).name));
Q_coefs = reshape(coefs_all, size(coefs_all, 1), size(coefs_all, 2), [], 3);

%combine Q and IB
IBcoefs_flat = reshape(IB_coefs, [], size(IB_coefs, 3), size(IB_coefs, 4));
Qcoefs_flat = reshape(Q_coefs, [], size(Q_coefs, 3), size(Q_coefs, 4));
combined_coefs = [IBcoefs_flat; Qcoefs_flat];


%%
paths = pathsetup('matchingsim');
% Load raw data
opts = struct;
opts.version = '121021';

[res, opts] = load_and_run(prob, opts);

[idxQ, idxIB] = reshapeidx(res.idx, res);

%% split into sub-classes and plot
coefs_means = {};
coefs_stds = {};
for i = 1:5
    coef_class = combined_coefs(res.idx == i, :, :);
    coefs_means{i} = squeeze(mean(coef_class, 1));
    coefs_stds{i} = squeeze(std(coef_class, [], 1));
end

savepath = fullfile(paths.figpath, 'history');

% plot
for i = 1:5
    figure('Position', [734,262,653,346]);
    imagesc(coefs_means{i}, 'YData', -1:-1:-5)
    caxis([-2, 2])
    colormap('redblue');
    colorbar('Position', [0.92,0.4,0.03,0.4], 'FontName', 'helvetica',...
        'FontSize', 14, 'FontAngle', 'italic')
    axis xy
    
    mymakeaxis('x_label', 'Factor', 'y_label', 'Past trials', 'xticks', 1:3, 'yticks', -1:-1:-5,...
        'xticklabels', {'Choice', 'Choice x Reward', 'Reward'})
    
    savename = fullfile(savepath, sprintf('history_coefs_split_by_class_class%d_122921.pdf', i));
    if ~exist(savename, 'file')
        saveas(gcf, savename);
        fprintf('File saved!\n')
    else
        fprintf('File exists, skipping save...\n');
    end
    
end



%% Plot the choice coefs
figure;
for i = 1:5
    subplot(2,3,i)
    imagesc(squeeze(coefs_reshaped(:,:,i,2)));
%     caxis([-1,1])
    colorbar
    
end