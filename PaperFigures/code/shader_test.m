%% Load the data
prob = 0;
[res, opts] = load_and_run(prob, 1);
load('svm_small_prob1_100221.mat')
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_092921_withMdl.mat')
% load(sprintf('%s/svmresults_from_pickle_092221_prob%.2f.mat', opts.rootdir, 1-opts.prob));
load(sprintf('%s/svmresults_from_pickle_100221_prob%.2f.mat', opts.rootdir, 1-opts.prob));

%

[idxQ, idxIB] = reshapeidx(res.idx, res);


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

Qoffsetall(Qoffsetall < -20) = 3;
IBoffsetall(IBoffsetall < -20) = 3;



features = [IBeffall IBlapseall IBslopeall IBoffsetall;
    Qeffall Qlapseall Qslopeall Qoffsetall];
featcopy = res.features;
% featcopy(:,3) = -featcopy(:,3);



features_norm1 = (features - mean(features, 1)) ./ std(features, [], 1);

V = pca(features_norm1);

labels = [idxIBall; idxQall];


figure;
for i = 1:4
    for j = 1:4
        subplot(4,4,4*(i-1) + j)
        hold on
        for k = 1:5
            plot(features(labels == k, i), features(labels == k, j), '.')
            plot(featcopy(res.idx == k, i), featcopy(res.idx == k, j), 'rx')
            
        end
    end
end

pred = Mdl.predict(features);
% pred = Mdls1{1}.predict(features);

% PC space
features_norm = (features - nanmean(res.features, 1)) ./ (nanstd(res.features, [], 1));
Xproj = features_norm1 * V;
figure;
hold on
for i = 1:5
    plot(Xproj(pred == i,1), Xproj(pred == i, 2), '.');
    
end



%%


Y = res.Y;
YminX = min(Y(:,1));
YmaxX = max(Y(:,1));
YminY = min(Y(:,2));
YmaxY = max(Y(:,2));

xcoords = linspace(-4, 4, 100);
ycoords = linspace(-4, 4, 100);
% [xx, yy] = meshgrid(xcoords, ycoords, ycoords);

domains = zeros(numel(xcoords), numel(ycoords));

% for each (x,y), find the Y entry that is closest
for i = 1:numel(xcoords)
    for j = 1:numel(ycoords)
        xval = xcoords(i);
        yval = ycoords(j);
        
        D = [xval yval] - res.Y;
        
        dist = D(:,1).^2 + D(:,2).^2;
        
%         Mdls1{1}.predict([xval yval]
        
        id = argmin(dist);
        
        domains(i,j) = res.idx(id);
        
        
    end
end

%% Plotting
figure;
imagesc(domains, 'XData', [min(xcoords) max(xcoords)], 'YData', [min(ycoords) max(ycoords)])
hold on

for i = 1:5
    plot(Xproj(idx == i, 2), Xproj(idx == i, 1), '.');
end

%%
figure;
copyproj = ftcopy * res.V;
hold on
for i = 1:5
    plot(copyproj(res.idx == i, 2), copyproj(res.idx == i, 1), '.');
end
%% With SVM predictions
xcoords = linspace(-4, 4, 100);
ycoords = linspace(-4, 4, 100);
ftcopy = res.features;
ftcopy(:,3) = -ftcopy(:,3);
minvals = min(ftcopy, [], 1);
maxvals = max(ftcopy, [], 1);
% grid = ngrid(linspace(minvals(1), maxvals(1), 10),...
%     linspace(minvals(2), maxvals(2), 10),...
%     linspace(minvals(3), maxvals(3), 10),...
%     linspace(minvals(4), maxvals(4), 10));
N = 5000;
x1 = unifrnd(minvals(1), maxvals(1), [N 1]);
x2 = unifrnd(minvals(2), maxvals(2), [N 1]);
x3 = unifrnd(minvals(3), maxvals(3), [N 1]);
x4 = unifrnd(minvals(4), 0, [N 1]);

X = [x1 x2 x3 x4];
idx = Mdls1{1}.predict(X);

Xnorm = (X - nanmean(X, 1)) ./ nanstd(X, [], 1);
Vnew = pca(Xnorm);

% idx2 = Mdls1{1}.predict(res.features_norm);

%project X to PC space
Xproj = Xnorm * Vnew;

for i = 1:numel(xcoords)
    for j = 1:numel(ycoords)
        xval = xcoords(i);
        yval = ycoords(j);
        
        D = [xval yval] - Xproj(:, 1:2);
        
        dist = D(:,1).^2 + D(:,2).^2;
        
        id = argmin(dist);
        
        domains(i,j) = idx(id);
        
        
    end
end



%%
% [xx, yy] = meshgrid(0:0.1:0.5, 0.2:0.3:1.4);