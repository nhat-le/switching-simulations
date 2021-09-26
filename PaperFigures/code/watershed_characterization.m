folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';
prob = 0;
switch prob
    case 0
        filename = 'opts_prob0.0-2021-09-25 20.52.mat';
    case 0.1
        filename = 'opts_prob0.1-2021-09-25 21.44.mat';
    case 0.2
        filename = 'opts_prob0.2-2021-09-25 21.57.mat';
    case 0.3
        filename = 'opts_prob0.3-2021-09-25 22.29.mat';
end

load(fullfile(folder, filename));
opts.save = 0;

%%
[idx, out] = run_watershed(opts);

% Rotation for idx
switch prob
    case 0
        idx = rotate(idx, [3, 2, 4]);
    case 0.1
        idx = rotate(idx, [6, 1, 3]);
        idx = rotate(idx, [4 5]);
    case 0.2
        idx = rotate(idx, [2, 1, 4]);
        idx = rotate(idx, [5, 3]);
    case 0.3
        idx = rotate(idx, [4, 1]);
        idx = rotate(idx, [2, 5, 3]);
end
opts.save = 0;

%% Regime characterization
features = out.features;
features1 = features(idx == 1,:);
features2 = features(idx == 2,:);
features3 = features(idx == 3,:);
features4 = features(idx == 4,:);
% features5 = features(idx == 5,:);


means = [mean(features1, 1) mean(features2, 1) mean(features3, 1), ...
    mean(features4, 1)];
stds = [std(features1, [], 1), std(features2, [], 1), std(features3, [], 1), ...
    std(features4, [], 1)];

make_feature_plot(1, means, stds, 'Efficiency', nan);
make_feature_plot(2, means, stds, 'Lapse', [0, 0.3]);
make_feature_plot(3, means, stds, 'Slope', [0, 20]);
make_feature_plot(4, -means, stds, 'Offset', nan);


%% Decoding
% Load simult sim dataset
load(sprintf('%s/svmresults_from_pickle_092221_prob%.2f.mat', opts.rootdir, 1-opts.prob));
[idxQ, idxIB] = reshapeidx(idx, out);

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

Qoffsetall(Qoffsetall < -20) = -3;
IBoffsetall(IBoffsetall < -20) = -3;



features = [IBeffall IBlapseall IBslopeall IBoffsetall;
    Qeffall Qlapseall Qslopeall Qoffsetall];
features_norm = (features - mean(features, 1)) ./ std(features, [], 1);

labels = [idxIBall; idxQall];

%shuffle
order = randperm(numel(labels));
labels_shuffled = labels(order);
features_shuffled = features(order,:);



%80% training, 20% testing
rng('shuffle');
ntrain = floor(numel(labels) * 0.8);
Xtrain = features_shuffled(1:ntrain,:);
ytrain = labels_shuffled(1:ntrain);
Xtest = features_shuffled(ntrain + 1:end,:);
ytest = labels_shuffled(ntrain + 1:end);

% t = templateLinear();
t = templateSVM('Standardize',true, 'KernelFunction', 'rbf');
Mdl = fitcecoc(Xtrain',ytrain,'Learners',t,'ObservationsIn','columns');
ypred = Mdl.predict(Xtest);

% Performance
perf = sum(ypred == ytest) / numel(ytest);

% Make the confusion matrix
counts = nan(5,5);
for i = 1:5
    for j = 1:5
        counts(i,j) = sum(ypred == i & ytest == j);

    end
end

confusion = counts ./ sum(counts, 1);

%%
figure;
hold off
confusionchart(ytest, ypred,'RowSummary','row-normalized');
% mymakeaxis('x_label', 'Predicted Class');
set(gca,'FontSize', 16, 'FontName', 'helvetica')


%%
[idxQ, idxIB] = reshapeidx(idx, out);
opts.N = numel(unique(idx));

cmap = brewermap(opts.N, 'Blues');
cmap = [cmap; 0.2 0.2 0.2];
h = produce_heatmap(idxIB, out.prewlst, out.pswitchlst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4, 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
ibfig = gcf;


produce_heatmap(idxQ, out.epslst, out.gammalst, 'clim', [0.5 opts.N+1.5], 'legendname', 'Performance regime', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'cmap', cmap, 'limits', [0.5, opts.N + 0.5]);
qfig = gcf;

if opts.save
    currdate = datetime;
    currdate.Format = 'yyyy-MM-dd HH.mm';
    currdate = string(currdate);    
    
    %Save Q plot
    filename = sprintf('%s/perfRegimeQ_permuted_prob%.1f-%s.pdf', opts.savepath, 1-opts.prob, currdate);
    if ~exist(filename, 'file')
        saveas(qfig, filename);
    end

    % Save IB plot
    filename = sprintf('%s/perfRegimeIB_permuted_prob%.1f-%s.pdf', opts.savepath, 1-opts.prob, currdate);
    if ~exist(filename, 'file')
        saveas(ibfig, filename);
    end
    
    fprintf('Files saved!\n');
    close all

end



function res = rotate(arr, order)
res = arr;
for i = 1:numel(order) - 1
    res(arr == order(i)) = order(i+1);
end

res(arr == order(end)) = order(1);

end

function make_feature_plot(classid, means, stds, featurename, ylims)
figure('Position', [440,376,320,422]);
paperaesthetics;
N = numel(means) / 4;
colors = brewermap(10, 'Blues');
% coltouse = colors(8, :);
coltouse = [8, 81, 156]/ 255;
errorbar(1:N, means(classid:4:end), stds(classid:4:end), 'o', 'MarkerFaceColor', coltouse,...
    'LineWidth', 2, 'MarkerEdgeColor', coltouse, 'Color', coltouse);

if ~isnan(ylims)
    ylim(ylims)
end

mymakeaxis('x_label', 'Class', 'y_label', featurename, 'xticks', 1:5)
end