%% Regime demarcation (with k-means clustering)
% Form the feature vectors
prob = 1;
savefig = 0;
expdate = '092321';
rootdir = fullfile('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/', expdate);

rng(127)
N = 5;


% Load Q-learning simulation results
filedir = sprintf('%s/EGreedyQLearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', rootdir, 1-prob, prob, expdate);
load(filedir);
Qeff_flat = reshape(efflist, [], 1);
Qlapse_flat = reshape(LapseL, [], 1);
Qslope_flat = reshape(PLslopelist, [], 1);
Qoffset_flat = reshape(PLoffsetlist, [], 1);
[Qxdim, Qydim] = size(efflist);

% Load inference-based simulation results
filedir = sprintf('%s/EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', rootdir, 1-prob, prob, expdate);
load(filedir);
IBeff_flat = reshape(efflist, [], 1);
IBlapse_flat = reshape(LapseL, [], 1);
IBslope_flat = reshape(PLslopelist, [], 1);
IBoffset_flat = reshape(PLoffsetlist, [], 1);
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
if savefig
    filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs/perfRegimeIB_prob%.1f.pdf', 1-prob);
    if ~exist(filename, 'file')
        saveas(gcf, filename);
    end
end


% subplot(122)
produce_heatmap(idxQ, epslst, gammalst, 'clim', [0.5 N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, N + 0.5]);
% Save IB plot
if savefig
    filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs/perfRegimeQ_prob%.1f.pdf', 1-prob);
    if ~exist(filename, 'file')
        saveas(gcf, filename);
    end
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

