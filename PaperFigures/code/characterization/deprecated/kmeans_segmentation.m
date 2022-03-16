%% Build the options for the run

% clear all
paths = pathsetup('matchingsim');
% where the Python simulation results are stored (model-free/inf-based
% simulations)
opts = struct;
opts.rootdir = paths.simdatapath;
opts.expfolder = '010422';
opts.prob = 0.8;
opts.seed = 11; %for p = 0.8, seed 10 (seed 4, 1.13.21)
opts.kclust = 6;
opts.clipmode = 2;
opts.rotations = {[1]};

% script operations
opts.save = 0;
opts.plotfeatures = 0;
opts.show_trans_functions = 0;
opts.method = 'kmeansraw'; %'multidim-watershed or 'watershed', 'watershed-tsne', or 'dbscan', or 'gmm'


colors = brewermap(12, 'Paired');
opts.cmap = colors([1,5,7,11,9,3],:); %permute to align with MATLAB default..


%% form the feature vectors
rng(opts.seed)
% rng(8)

[out,opts] = load_data(opts);

all_transfuncs = [];
xvals = 1:15;

for i = 1:size(out.features, 1)
    mu = -out.features(i, 4);
    sigma = out.features(i, 3);
    lapse = out.features(i, 2);
    transfunc = mathfuncs.sigmoid(xvals, mu, sigma, lapse);
    all_transfuncs(end+1,:) = transfunc;
    
end

idx = kmeans(all_transfuncs, opts.kclust); %, 'Distance', 'cityblock');
idx = rotate_idx(idx, opts.rotations);
out.idx = idx;

figure;
for clustid = 1:max(idx)
    subplot(2, ceil(opts.kclust/2), clustid)
    hold on
    transfunc_sub = all_transfuncs(idx == clustid, :)';
    plot(transfunc_sub(:,:), 'k', 'LineWidth', 0.25);
    plot(mean(transfunc_sub, 2), 'r', 'Color', 'r', 'LineWidth', 2);
    
    ylim([0, 1])
    mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)', ...
        'yticks', 0:0.2:1, 'xticks', 0:5:15);
end
clustfig = gcf;


%
[idxQ, idxIB] = reshapeidx(idx, out);
opts.N = numel(unique(idx));

cmap = opts.cmap;
% cmap = brewermap(opts.N, 'Blues');
cmap = [cmap; 0.2 0.2 0.2];
h = produce_heatmap(idxIB, out.prewlst, out.pswitchlst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
ibfig = gcf;


produce_heatmap(idxQ, out.epslst, out.gammalst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
qfig = gcf;


%% Plot all the transition functions in one plot
figure;
currdate = datetime;
currdate.Format = 'yyyy-MM-dd HH.mm';
currdate = string(currdate);    
plot(all_transfuncs', 'k', 'LineWidth', 0.25);
mymakeaxis('x_label', 'Trials in block', 'y_label', 'P(Correct)', 'xticks', 0:5:15)
filename = sprintf('%s/all_transfuncs-%s-%s.pdf', opts.figsavepath, ...
    currdate, opts.method);
exportgraphics(gcf,filename,'ContentType','vector')

%% PC space projection
[~,~,V] = svd(out.features_norm);
out.V = V;
out.Y = out.features_norm * V;

colors = brewermap(7, 'Set1');
colors = colors([2,1,5,6,4,3,7],:);

figure()
    
hold on
for i = 1:numel(unique(idx))
    col2use = colors(i,:);
      
    plot(out.Y(idx == i, 1), out.Y(idx == i, 2), 'o', ...
        'MarkerFaceColor', col2use, 'MarkerEdgeColor', 'w', 'MarkerSize', 7)
end
mymakeaxis('x_label', 'PC1', 'y_label', 'PC2')


%% Validation for K selection
rng(opts.seed)
rng(15)
Kvals = 1:10;
eva = evalclusters(all_transfuncs,'kmeans','gap','KList',Kvals, ...
    'Distance', 'sqEuclidean');
figure;
plot(Kvals, eva.CriterionValues, 'LineWidth', 2)
xlim([1,10])
mymakeaxis('x_label', 'K', 'y_label', 'Silhouette metric', 'xticks', Kvals)
kvalfig = gcf;


%%
% res = [];
% for seed = 1:100
%     rng(seed)
%     Kvals = 1:10;
%     eva = evalclusters(all_transfuncs,'kmeans','silhouette','KList',Kvals);
%     res(seed) = eva.OptimalK;
%     figure;
%     plot(Kvals, eva.CriterionValues, 'LineWidth', 2)
%     xlim([1,10])
%     mymakeaxis('x_label', 'K', 'y_label', 'Silhouette metric', 'xticks', Kvals)
%     kvalfig = gcf;
% end

%% Dim-reduction
% Y = tsne(all_transfuncs, 'Perplexity', 100);
% figure;
% scatter(Y(:,1), Y(:,2), [], idx_clust)
% tsne on the raw data gave poor results for visualing and does not converge


%% Save if requested

if opts.save
    currdate = datetime;
    currdate.Format = 'yyyy-MM-dd HH.mm';
    currdate = string(currdate);    
    
    %Save Q plot
    filename = sprintf('%s/perfRegimeQ_prob%.1f-%s-%s.pdf', opts.figsavepath, ...
        1-opts.prob, currdate, opts.method);
    opts.fig = qfig;
    opts.type = 'pdf';
    ioutils.savesafe(filename, opts);
%     if ~exist(filename, 'file')
%         saveas(qfig, filename);
%     end

    % Save IB plot
    filename = sprintf('%s/perfRegimeIB_prob%.1f-%s-%s.pdf', opts.figsavepath,...
        1-opts.prob, currdate, opts.method);
    opts.fig = ibfig;
    ioutils.savesafe(filename, opts);
%     if ~exist(filename, 'file')
%         saveas(ibfig, filename);
%     end
    
    % Save cluster panel plot
    filename = sprintf('%s/clusters_prob%.1f-%s-%s.pdf', opts.figsavepath,...
        1-opts.prob, currdate, opts.method);
    opts.fig = clustfig;
    exportgraphics(clustfig,filename,'ContentType','vector')
%     ioutils.savesafe(filename, opts);
%     if ~exist(filename, 'file')
%         saveas(ibfig, filename);
%     end
    
    % Save IB plot
%     filename = sprintf('%s/kval_prob%.1f-%s-%s.pdf', opts.figsavepath,...
%         1-opts.prob, currdate, opts.method);
%     opts.fig = kvalfig;
%     ioutils.savesafe(filename, opts);
%     if ~exist(filename, 'file')
%         saveas(ibfig, filename);
%     end
%     
    % Save mds fig
%     filename = sprintf('%s/mds_prob%.1f-%s-%s.pdf', opts.figsavepath, ...
%         1-opts.prob, currdate, opts.method);
%     if ~exist(filename, 'file')
%         saveas(mdsfig, filename);
%     end
    
    % Save options
    filename = sprintf('%s/opts_prob%.1f-%s-%s.mat', opts.datasavepath,...
        1-opts.prob, currdate, opts.method);
    opts.fig = [];
    if ~exist(filename, 'file')
        save(filename, 'opts', 'out');
    else
        fprintf('File exists, skipping save: %s\n', filename)
    end
    
    fprintf('Files saved!\n');
    close all

end







