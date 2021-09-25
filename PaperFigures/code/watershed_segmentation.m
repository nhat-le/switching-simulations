% Build the options for the run
opts = struct;
opts.nbinhist = 25;
opts.imhmin = 1;
opts.kernelsize = 3;
opts.N = 6;
opts.prob = 1;

%% form the feature vectors
out = load_data(opts);

%
D = pdist(out.features_norm);
D = squareform(D);
[Y,e] = cmdscale(D,2);

figure;
plot(Y(:,1), Y(:,2), '.');


%%
f = figure;
h = histogram2(Y(:,1), Y(:,2), opts.nbinhist);
vals = h.Values;
binX = h.XBinEdges;
binY = h.YBinEdges;

% Lowpass filter
vals = conv2(double(vals), double(ones(opts.kernelsize, opts.kernelsize)), 'same');
vals = imhmin(-vals, opts.imhmin);

ordersX = [];
ordersY = [];
for i = 1:size(Y,1)
    ordersX(i) = find_order_in_arr(binX, Y(i,1)) - 1;
    ordersY(i) = find_order_in_arr(binY, Y(i,2)) - 1;
end

close(f)

% Watershed!
L = watershed(vals);
figure;
subplot(121)
imagesc(vals)
colormap gray
caxis([-200 100])
subplot(122)
vals(L == 0) = 100;
imagesc(vals)
caxis([-200 100])



%%
for i = 1:size(Y,1)
    labels(i) = L(ordersX(i), ordersY(i));
end

% Linterest = [2,5,4,7,8];
Linterest = [1 2 3 4 5];
idx = labels;
for i = 1:numel(Linterest)
    idx(labels == Linterest(i)) = i;
end
idx(labels == 0) = nan;

figure;
hold on
plot(Y(idx==1,1), Y(idx==1, 2), '.')
plot(Y(idx==2,1), Y(idx==2, 2), '.')
plot(Y(idx==3,1), Y(idx==3, 2), '.')
plot(Y(idx==4,1), Y(idx==4, 2), '.')
plot(Y(idx==5,1), Y(idx==5, 2), '.')
% plot(Y(idx>5,1), Y(idx>5, 2), '.')




%% Plot the embedded data
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

%% Post-processing
[idxQ, idxIB] = reshapeidx(idx, out);

cmap = brewermap(opts.N, 'Blues');
cmap = [cmap; 0.2 0.2 0.2];
h = produce_heatmap(idxIB, out.prewlst, out.pswitchlst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);

% Save Q plot
filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs/perfRegimeIB_prob%.1f.pdf', 1-opts.prob);
if ~exist(filename, 'file')
    saveas(gcf, filename);
end


% subplot(122)
produce_heatmap(idxQ, out.epslst, out.gammalst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
% Save IB plot
filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs/perfRegimeQ_prob%.1f.pdf', 1-opts.prob);
if ~exist(filename, 'file')
    saveas(gcf, filename);
end







