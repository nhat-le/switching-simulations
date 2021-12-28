paths = pathsetup('matchingsim');
expfitdate = '121821';
rootdir = fullfile(paths.blockhmmfitpath, expfitdate);
folders = dir(fullfile(rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
mdltypes = 1:2;
mdlids = 1:10;
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.python_assist = 0;
opts.effmethod = 'sim';

[~, aggmeans_native, aggparams_native] = load_params(folders, opts);

% concat and do k-means
aggparams_all = cell2mat(aggparams_native)';
aggmeans_all = cell2mat(aggmeans_native');

K = 5;
colors = brewermap(K, 'Set1');



%% visualize
figure;
for i = 1:3
    for j = i+1:4
        subplot(3,4, (i-1)*4+j);
        plot(aggparams_all(i,:), aggparams_all(j,:), '.')
    
    end
end


%% kmeans (plotting raw)
aggparams_all(aggparams_all > 20) = 20;

for i = 15
    rng(i)
    K = 4;
    aggparams_norm = (aggparams_all - mean(aggparams_all)) ./ std(aggparams_all, [], 1);
    idx = kmeans(aggparams_norm(:,1:3), K);

    figure;
    for j = 1:K
        subplot(2,3,j)
        aggmeans_sub = aggmeans_all(idx == j,:);
        plot(aggmeans_sub', 'b');


    end
end



%% kmeans (plotting smoothed)
aggparams_all(aggparams_all > 20) = 20;
for i = 15
    rng(i)
    K = 4;
    aggparams_norm = (aggparams_all - mean(aggparams_all)) ./ std(aggparams_all, [], 1);
    idx = kmeans(aggparams_norm(:,1:3), K);

    figure;
    hold on
    for j = 1:K
        subplot(2,3,j)
        hold on
        param_sub = aggparams_all(idx == j,:);
        for k = 1:size(param_sub, 1)
            param = param_sub(k,:);
            yvals = mathfuncs.sigmoid(1:20, param(1), param(2), param(3));
            plot(yvals, 'b');
        end


    end
end


%% with visualization
rng(15)
K = 5;
aggparams_norm = (aggparams_all - mean(aggparams_all)) ./ std(aggparams_all, [], 1);
% idx = kmeans(aggparams_norm(:,1:3), K);
labels = {'offset', 'slope', 'lapse', 'eff'};


figure;
colors = brewermap(K, 'Set1');

for xcol = 1:4
    for ycol = 1:4
        for k = 1:K
            subparams = aggparams_all(idx == k, :);
            subparams(:,1) = log10(subparams(:,1));
            subplot(4,4, (xcol - 1) * 4 + ycol);
            plot(subparams(:,xcol), subparams(:,ycol),'.', 'Color', colors(k,:), 'MarkerSize', 10);
            xlabel(labels{xcol})
            ylabel(labels{ycol})
            hold on
            
        end
        
    end
end


%% project to pc space
aggcopy = aggparams_all;
aggcopy(:,1) = log10(aggcopy(:,1));

V = pca(aggcopy);
Yproj = aggcopy * V;







%% plotting raw
rng(15)
K = 4;
aggparams_norm = (aggparams_all - mean(aggparams_all)) ./ std(aggparams_all, [], 1);
% idx = kmeans(aggparams_norm(:, 1:3), K);
labels = {'offset', 'slope', 'lapse', 'eff'};

figure;
for j = 1:K
    subplot(2,3,j)
    aggmeans_sub = aggmeans_all(idx == j,:);
    plot(aggmeans_sub', 'Color', colors(j,:));
    hold on
    plot(mean(aggmeans_sub, 1), 'b', 'Color', 'k', 'LineWidth', 2);

end


%% kmeans (plotting smoothed)
aggparams_all(aggparams_all > 20) = 20;
for i = 15
    rng(i)
    K = 4;
    aggparams_norm = (aggparams_all - mean(aggparams_all)) ./ std(aggparams_all, [], 1);
%     idx = kmeans(aggparams_norm(:,1:3), K);

    figure;
    hold on
    for j = 1:K
        subplot(2,3,j)
        hold on
        param_sub = aggparams_all(idx == j,:);
        for k = 1:size(param_sub, 1)
            param = param_sub(k,:);
            yvals = mathfuncs.sigmoid(1:20, param(1), param(2), param(3));
            plot(yvals, 'Color', colors(j, :));
        end


    end
end


%% manual clustering based on features
slopes = aggparams_all(:,2);
offsets = aggparams_all(:,1);
lapses = aggparams_all(:,3);
clust4_idx = (offsets < 2 & lapses < 0.2); %purple
clust3_idx = (offsets > 4.8); %green
clust2_idx = (lapses > 0.2 & offsets < 5); %blue
% clust1_idx = 1 - clust2_idx - clust3_idx - clust4_idx; %red
clust1_idx = (~clust4_idx & ~clust3_idx & ~clust2_idx & slopes < 3 & offsets > 1);
clust5_idx = 1 - clust1_idx - clust2_idx - clust3_idx - clust4_idx;


clust_all_idx = {clust1_idx, clust2_idx, clust3_idx, clust4_idx};
idx = clust1_idx + clust2_idx * 2 + clust3_idx * 3 + clust4_idx * 4 + clust5_idx * 5;

% plotting raw
figure;
for j = 1:K
    subplot(2,3,j)
    aggmeans_sub = aggmeans_all(idx == j,:);
    plot(aggmeans_sub', 'Color', colors(j,:));
    hold on
    plot(mean(aggmeans_sub, 1), 'b', 'Color', 'k', 'LineWidth', 2);

end
% for j = 1:4
%     idx_clust = clust_all_idx{j};
%     subplot(2,2,j)
%     hold on
%     param_sub = aggparams_all(logical(idx_clust),:);
%     for k = 1:size(param_sub, 1)
%         param = param_sub(k,:);
%         yvals = mathfuncs.sigmoid(1:20, param(1), param(2), param(3));
%         plot(yvals, 'Color', colors(j, :));
%     end
% end

%%
figure;
for j = 1:4
    idx_clust = clust_all_idx{j};
    subplot(2,3,j)
    aggmeans_sub = aggmeans_all(logical(idx_clust),:);
    plot(aggmeans_sub', 'b');


end







    
    