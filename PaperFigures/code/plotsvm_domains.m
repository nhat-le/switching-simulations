%% Load the data
prob = 0;
[res, opts] = load_and_run(prob, 1);
load('svm_small_prob1_100221.mat')


% svm predictions
xcoords = linspace(-4, 4, 100);
ycoords = linspace(-4, 4, 100);
ftcopy = res.features;
ftcopy(:,3) = -ftcopy(:,3);
minvals = min(ftcopy, [], 1);
maxvals = max(ftcopy, [], 1);

N = 5000;
x1 = unifrnd(minvals(1), maxvals(1), [N 1]);
x2 = unifrnd(minvals(2), maxvals(2), [N 1]);
x3 = unifrnd(minvals(3), maxvals(3), [N 1]);
x4 = unifrnd(minvals(4), 0, [N 1]);

X = [x1 x2 x3 x4];
idx = Mdl.predict(X);

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