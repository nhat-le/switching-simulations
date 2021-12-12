%% Build the options for the run

clear all

% where the Python simulation results are stored (model-free/inf-based
% simulations)
opts = struct;
opts.rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata';
opts.expfolder = '121021';
opts.nbinhist = 30;
opts.imhmin = 3;
opts.kernelsize = 3;
opts.prob = 1;
opts.save = 1;
opts.seed = 3;
opts.rotations = {[3, 5], [2, 4]};
opts.plotfeatures = 0;


rng(opts.seed);

%% form the feature vectors
[out,opts] = load_data(opts);

D = pdist(out.features_norm);
D = squareform(D);
[Y,e] = cmdscale(D,2);

V = pca(out.features_norm);
Xproj = out.features_norm * V;

% figure;
% plot(Y(:,1), Y(:,2), '.');
% figure;
% plot(Xproj(:,1), Xproj(:,2), '.');



%%
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


%%
% labels = assign_labels(Y, L, binX, binY, 'Lfill', Lfill);
labels = assign_labels(Y, L, binX, binY, 'nearest');
%
nLlabels = unique(Lfill);
fprintf('unique class labels: %d\n', numel(nLlabels));
idx = labels;

figure;
hold on
for i = nLlabels'
    plot(Y(idx == i,1), Y(idx == i, 2), '.')
end

mdsfig = gcf;


%% Plot the embedded data
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
[idxQ, idxIB] = reshapeidx(idx, out);
opts.N = numel(nLlabels);

cmap = brewermap(opts.N, 'Blues');
cmap = [cmap; 0.2 0.2 0.2];
h = produce_heatmap(idxIB, out.prewlst, out.pswitchlst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
ibfig = gcf;


produce_heatmap(idxQ, out.epslst, out.gammalst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
qfig = gcf;



%% Save if requested

if opts.save
    currdate = datetime;
    currdate.Format = 'yyyy-MM-dd HH.mm';
    currdate = string(currdate);    
    
    %Save Q plot
    filename = sprintf('%s/perfRegimeQ_prob%.1f-%s.pdf', opts.figsavepath, 1-opts.prob, currdate);
    if ~exist(filename, 'file')
        saveas(qfig, filename);
    end

    % Save IB plot
    filename = sprintf('%s/perfRegimeIB_prob%.1f-%s.pdf', opts.figsavepath, 1-opts.prob, currdate);
    if ~exist(filename, 'file')
        saveas(ibfig, filename);
    end
    
    % Save mds fig
    filename = sprintf('%s/mds_prob%.1f-%s.pdf', opts.figsavepath, 1-opts.prob, currdate);
    if ~exist(filename, 'file')
        saveas(mdsfig, filename);
    end
    
    % Save options
    filename = sprintf('%s/opts_prob%.1f-%s.mat', opts.datasavepath, 1-opts.prob, currdate);
    if ~exist(filename, 'file')
        save(filename, 'opts');
    end
    
    fprintf('Files saved!\n');
    close all

end







