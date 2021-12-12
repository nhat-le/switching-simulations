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


%% Decoding results for 20-block sessions
% Load simult sim dataset
svmrootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/svm/';
svmdate = '092221';
prob = 1;
load(sprintf('%s/%s/svmresults_from_pickle_%s_prob%.2f.mat', svmrootdir, svmdate, svmdate, 1-prob));
% idxQresize = ceil(imresize(idxQ, [11 11]));
% idxIBresize = floor(imresize(idxIB, [11 11]));

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


% Qoffsetall(Qoffsetall < -20) = -20;
% IBoffsetall(IBoffsetall < -20) = -20;


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

% save('svm_Mdl.mat', 'Mdl', 'idxQ' , 'idxIB');

%% Visualize clusters and labels
figure;
hold on
for i = 1:4
    for j = 1:4
        subplot(4,4,(i-1)*4+j)
        hold on
        plot(features(labels==1,i), features(labels==1,j), '.'); %blue
        plot(features(labels==2,i), features(labels==2,j), '.'); %red
        plot(features(labels==3,i), features(labels==3,j), '.'); %yellow
        plot(features(labels==4,i), features(labels==4,j), '.'); %purple
        plot(features(labels==5,i), features(labels==5,j), '.'); %green

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


%%
figure;
hold off
confusionchart(ytest, ypred,'RowSummary','row-normalized');
% mymakeaxis('x_label', 'Predicted Class');
set(gca,'FontSize', 16, 'FontName', 'helvetica')









function visualize_kmeans(i, j, idx, features_norm)
subplot(4,4,(i-1) * 4 + j)
plot(features_norm(idx == 1, i), features_norm(idx == 1, j), '.');
hold on
plot(features_norm(idx == 2, i), features_norm(idx == 2, j), '.');
plot(features_norm(idx == 3, i), features_norm(idx == 3, j), '.');

end

function make_feature_plot(classid, means, stds, featurename, ylims)
figure('Position', [440,376,320,422]);
paperaesthetics;
colors = brewermap(10, 'Blues');
% coltouse = colors(8, :);
coltouse = [8, 81, 156]/ 255;
errorbar(1:5, means(classid:4:end), stds(classid:4:end), 'o', 'MarkerFaceColor', coltouse,...
    'LineWidth', 2, 'MarkerEdgeColor', coltouse, 'Color', coltouse);

if ~isnan(ylims)
    ylim(ylims)
end

mymakeaxis('x_label', 'Class', 'y_label', featurename, 'xticks', 1:5)
end


function res = rotate(arr, order)
res = arr;
for i = 1:numel(order) - 1
    res(arr == order(i)) = order(i+1);
end

res(arr == order(end)) = order(1);

end

