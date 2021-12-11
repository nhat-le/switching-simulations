load('simdata/EGreedyqlearningAgent-withCorr-doublesigmoid-prob0.00to1.00-100221.mat');

opts.prob = 1;
opts.filestem{1} = 'EGreedyqlearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-100221.mat';
opts.filestem{2} = 'EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-100221.mat';
[out,opts] = load_data(opts);


%%
opts2.prob = 1;
[out2,opts2] = load_data(opts2);

%eff_flat = reshape(efflist



%%

[res1, out1] = load_and_run(0);
[idxQ, idxIB] = reshapeidx(res1.idx, res1);

%% 
featnorm1 = (out.features - mean(out2.features, 1)) ./ std(out2.features, [], 1);
featnorm2 = (out2.features - mean(out2.features, 1)) ./ std(out2.features, [], 1);

Y1 = featnorm1 * res1.V;
Y2 = featnorm2 * res1.V;

figure;
subplot(121)
plot(Y1(:,1), Y1(:,2), '.')

subplot(122)
plot(Y2(:,1), Y2(:,2), '.')

[idxQsim, idxIBsim] = reshapeidx(1:650, res1);


%% svm on small data
labels = res1.idx';
order = randperm(numel(labels));

% order = randsample(1:numel(labels), numel(labels), true);
labels_shuffled = labels(order);
features_shuffled = res1.features(order,:);

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
N = numel(unique(ypred));
fprintf('Number of clusters = %d\n', N);
counts = nan(N,N);
for i = 1:N
    for j = 1:N
        counts(i,j) = sum(ypred == i & ytest == j);

    end
end

% Draw the decision boundaries









