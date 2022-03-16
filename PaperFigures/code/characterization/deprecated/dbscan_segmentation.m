%% Build the options for the run

% clear all
paths = pathsetup('matchingsim');
% where the Python simulation results are stored (model-free/inf-based
% simulations)
opts = struct;
opts.rootdir = paths.simdatapath;
opts.expfolder = '010422';
opts.nbinhist = 50;
opts.imhmin = 4;
opts.kernelsize = 10;
opts.prob = 1;
opts.seed = 130;
opts.rotations = nan; %{[4, 2, 1], [5, 3]};

% script operations
opts.save = 0;
opts.plotfeatures = 0;
opts.show_trans_functions = 0;

% parameters for tsne
opts.perplexity = 50;

% parameters for dbscan
opts.minpts = 40;
opts.epsilon = 0.4;


rng(opts.seed);

opts.method = 'dbscan'; %'multidim-watershed or 'watershed', 'watershed-tsne', or 'dbscan', or 'gmm'

% [out,opts] = load_data(opts);

% Y = tsne(out.features_norm);


%% form the feature vectors
[out,opts] = load_data(opts);

D = pdist(out.features_norm);
D = squareform(D);


if strcmp(opts.method, 'watershed-tsne')
    [Ypc,e] = cmdscale(D,2);
    Y = tsne(out.features_norm, 'Perplexity', opts.perplexity, 'Distance', 'cosine');
else
    [Y,e] = cmdscale(D,2);
end

V = pca(out.features_norm);
Xproj = out.features_norm * V;


%% multi-dimensional watershed
switch opts.method
    case 'multidim-watershed'
        nbin = 10;
        [counts, edges, mid, loc] = histcn(out.features_norm, nbin, nbin, nbin, nbin);
        filtsize = 2;
        filter = ones(filtsize,filtsize,filtsize,filtsize) / filtsize^4;%(counts * 0 + 1) / numel(counts);
        counts2 = convn(counts, filter, 'same');
        L = watershed(counts2);

        %TODO: assign cluster id to each training point


    case 'gmm'
        GMMModel = fitgmdist(out.features_norm, 4, 'Replicates', 10);
        idx = cluster(GMMModel, out.features_norm);
        labels_name = unique(idx)';



    case {'watershed', 'watershed-tsne'}
        f = figure;
        h = histogram2(Y(:,1), Y(:,2), opts.nbinhist);
        vals = h.Values;
        binX = h.XBinEdges;
        binY = h.YBinEdges;

        % Lowpass filter
        vals = conv2(double(vals), double(ones(opts.kernelsize, opts.kernelsize)), 'same');
        vals = imhmin(-vals, opts.imhmin);


        % vals = imhmax(vals, 30);


        close(f)

        % Watershed!
        L = watershed(vals);
        Lfill = remove_borders(L);

        figure;
        subplot(131)
        imagesc(vals')
        colormap gray
        caxis([-200 100])
        axis xy

        subplot(132)
        valsold = vals;
        vals(L == 0) = 100;
        imagesc(vals')
        caxis([-200 100])
        axis xy


        subplot(133)
        imagesc(Lfill', 'AlphaData', -valsold' / -min(valsold(:)) * 3);
        colormap jet
        axis xy

        idx = assign_labels(Y, L, binX, binY, 'nearest');

        labels_name = unique(idx);
        fprintf('Num of unique classes = %d\n', numel(unique(idx)));


    case 'dbscan'
        idx = dbscan(out.features_norm, opts.epsilon, opts.minpts);
        fprintf('unique class labels: %d\n', numel(unique(idx)));
        labels_name = unique(idx)';
        
        % since we have a -1 label, shift it to the max label instead
        idxcopy = idx;
        idx(idxcopy == -1) = max(labels_name) + 1;
        labels_name = unique(idx)';

end

% labels = assign_labels(Y, L, binX, binY, 'Lfill', Lfill);
%
% nLlabels = unique(Lfill);
% idx = labels;

%minpts = 10, eps = 0.5 for p = 1
%minpts = 10, eps= 0.4 for p = 0.9
%minpts = 10, eps= 0.3 for p = 0.8
%minpts = 15, eps= 0.35 for p = 0.8
%minpts = 20, eps= 0.4 for p = 0.8
%minpts = 20, eps = 0.42 for p = 0.7
%minpts = 15, eps = 0.4 for p = 0.7
%minpts = 10, eps = 0.35 for p = 0.7


%1.4.22 dataset



% kD = pdist2(out.features_norm, out.features_norm, 'euc', 'Smallest', minpts);
% plot(sort(kD(end,:)));

%%
% figure;
% hold on;
% nLabels = numel(labels_name);
% for i = 1:numel(labels_name) %nLlabels'
%     idplot = labels_name(i);
%     subplot(nLabels,2,i)
%     plot(Ypc(idx == idplot,1), Ypc(idx == idplot, 2), '.')
%     
%     subplot(nLabels,2, i+nLabels)
%     plot(Y(idx == idplot,1), Y(idx == idplot, 2), '.')
% end
% 
% mdsfig = gcf;


%%
figure;
hold on;
l = [];
names = {};
for i = labels_name %nLlabels'
%     plot(Ypc(idx == i,1), Ypc(idx == i, 2), '.')
    
    hold on
    l(end+1) = plot(Y(idx == i,1), Y(idx == i, 2), '.');
    names{end+1} = num2str(i);
end

mdsfig = gcf;
legend(l, names);


%% Plot the embedded data
% idxcopy = idx;
% idx(idxcopy == -1) = 5;

if opts.plotfeatures
    figure;
    hold on
    for i = 1:4
        for j = 1:4
            subplot(4,4,(i-1)*4+j)
            hold on
            plot(out.features(idx==1,i), out.features(idx==1,j), '.'); %blue
            plot(out.features(idx==2,i), out.features(idx==2,j), '.'); %red
            plot(out.features(idx==3,i), out.features(idx==3,j), '.'); %yellow
            plot(out.features(idx==4,i), out.features(idx==4,j), '.'); %purple
            plot(out.features(idx==5,i), out.features(idx==5,j), '.'); %green

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
end

%% Post-processing
% [idxQ, idxIB] = reshapeidx(idx, out);
% opts.N = numel(unique(idx));

[idxQ, idxIB] = reshapeidx(idx_clust, out);
opts.N = numel(unique(idx_clust));

cmap = brewermap(opts.N, 'Blues');
cmap = [cmap; 0.2 0.2 0.2];
h = produce_heatmap(idxIB, out.prewlst, out.pswitchlst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
ibfig = gcf;


produce_heatmap(idxQ, out.epslst, out.gammalst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
qfig = gcf;


%% Clustering using the transition functions
all_transfuncs = [];
xvals = 1:15;

for i = 1:size(out.features, 1)
    mu = -out.features(i, 4);
    sigma = out.features(i, 3);
    lapse = out.features(i, 2);
    transfunc = mathfuncs.sigmoid(xvals, mu, sigma, lapse);
    all_transfuncs(end+1,:) = transfunc;
    
end

figure;
idx_clust = kmeans(all_transfuncs, 6, 'Distance', 'cityblock');
for clustid = 1:max(idx_clust)
    subplot(2,4,clustid)
    hold on
    plot(all_transfuncs(idx_clust == clustid, :)', 'k');
    
    ylim([0, 1])
end

%%
rng(1)
Kvals = 1:15;
eva = evalclusters(all_transfuncs,'kmeans','silhouette','KList',Kvals);
figure;
plot(Kvals, eva.CriterionValues)

%% Visualize the clustered transition functions
figure;
for clustid = 1:numel(unique(idx))
    subplot(2,4,clustid)
    hold on
    features_clust = out.features(idx == clustid, :);
    xvals = 1:15;
    for i = 1:size(features_clust, 1)
        transfunc = mathfuncs.sigmoid(xvals, -features_clust(i,4),...
            features_clust(i,3), features_clust(i,2));
        plot(transfunc, 'k')
    end
    
    ylim([0, 1])
end


%% Save if requested

if opts.save
    currdate = datetime;
    currdate.Format = 'yyyy-MM-dd HH.mm';
    currdate = string(currdate);    
    
    %Save Q plot
    filename = sprintf('%s/perfRegimeQ_prob%.1f-%s-%s.pdf', opts.figsavepath, ...
        1-opts.prob, currdate, opts.method);
    if ~exist(filename, 'file')
        saveas(qfig, filename);
    end

    % Save IB plot
    filename = sprintf('%s/perfRegimeIB_prob%.1f-%s-%s.pdf', opts.figsavepath,...
        1-opts.prob, currdate, opts.method);
    if ~exist(filename, 'file')
        saveas(ibfig, filename);
    end
    
    % Save mds fig
    filename = sprintf('%s/mds_prob%.1f-%s-%s.pdf', opts.figsavepath, ...
        1-opts.prob, currdate, opts.method);
    if ~exist(filename, 'file')
        saveas(mdsfig, filename);
    end
    
    % Save options
    filename = sprintf('%s/opts_prob%.1f-%s-%s.mat', opts.datasavepath,...
        1-opts.prob, currdate, opts.method);
    if ~exist(filename, 'file')
        save(filename, 'opts');
    end
    
    fprintf('Files saved!\n');
    close all

end







