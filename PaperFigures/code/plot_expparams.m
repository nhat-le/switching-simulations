% load expdata/expfit_params.mat
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params.mat');
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_092921.mat');

Nsess = 20; %plot only this many sessions
bluecol = [8,81,156]/255;

%% Eff
figure;
% plot(expeff_all', 'Color', [55,126,184]/255)
effmeans = nanmean(expeff_all, 1);
efferr = nanstd(expeff_all, [], 1) / sqrt(size(expeff_all, 1));

hold on
% plot(nanmean(expeff_all, 1), 'Color', [228,26,28]/255, 'LineWidth', 2)
errorbar(1:Nsess, effmeans(1:Nsess), efferr(1:Nsess), 'o', ...,
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
mymakeaxis('x_label', 'Session', 'y_label', 'Performance', 'xytitle', 'N = 15',...
    'yticks', 0.5:0.1:0.7)

% Fit the eff mean
xvals = 1:Nsess;
yvals = effmeans(1:Nsess);
p = [1, 10, 0.2];
pfit = fminsearch(@(p) losssigmoid(p, xvals, yvals), p);
ypred = 0.5 + (1 - 0.5 - pfit(3)) ./ (1 + exp(-pfit(1) * (xvals - pfit(2))));
ypred2 = 0.5 + (1 - 0.5 - p(3)) ./ (1 + exp(-p(1) * (xvals - p(2))));

plot(xvals, ypred, 'LineWidth', 2, 'Color', bluecol);


%% offset
figure;
expoffsets_all(expoffsets_all > 100) = nan;
% plot(expeff_all', 'Color', [55,126,184]/255)
offsetmeans = nanmean(expoffsets_all, 1);
offseterr = nanstd(expoffsets_all, [], 1) / sqrt(size(expoffsets_all, 1));

% Fit the offset mean
xvals = 1:Nsess;
yvals = offsetmeans(1:Nsess);
p = [1, 1, 2];
lossfun(p, xvals, yvals)
pfit = fminsearch(@(p) lossfun(p, xvals, yvals), p);
ypred = pfit(1) * exp(-pfit(2)*xvals) + pfit(3);


% plot(nanmean(expeff_all, 1), 'Color', [228,26,28]/255, 'LineWidth', 2)
errorbar(1:Nsess, offsetmeans(1:Nsess), offseterr(1:Nsess), 'o', ...,
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
hold on
plot(xvals, ypred, 'LineWidth', 2, 'Color', bluecol)
mymakeaxis('x_label', 'Session', 'y_label', 'Offset', 'xytitle', 'N = 15')


%% Slopes
figure;
% plot(expeff_all', 'Color', [55,126,184]/255)
slopemeans = nanmean(expslopes_all, 1);
slopeerr = nanstd(expslopes_all, [], 1) / sqrt(size(expslopes_all, 1));

hold on
% plot(nanmean(expeff_all, 1), 'Color', [228,26,28]/255, 'LineWidth', 2)
errorbar(1:Nsess, slopemeans(1:Nsess), slopeerr(1:Nsess), 'o', ...,
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
mymakeaxis('x_label', 'Session', 'y_label', 'Slope', 'xytitle', 'N = 15')

% Linear regression 
xvals = 1:Nsess;
yvals = slopemeans(1:Nsess);
X = [ones(length(xvals),1) xvals'];
b = X\yvals';
ypred = b(1) + b(2) * xvals;
plot(xvals, ypred, 'LineWidth', 2, 'Color', bluecol);


%% lapses
figure;
% plot(expeff_all', 'Color', [55,126,184]/255)
lapsemean = nanmean(explapses_all, 1);
lapseerr = nanstd(explapses_all, [], 1) / sqrt(size(explapses_all, 1));

% Fit the lapse mean
xvals = 1:Nsess;
yvals = lapsemean(1:Nsess);
p = [1, 1, 0.1];
lossfun(p, xvals, yvals)
pfit = fminsearch(@(p) lossfun(p, xvals, yvals), p);
ypred = pfit(1) * exp(-pfit(2)*xvals) + pfit(3);



hold on
% plot(nanmean(expeff_all, 1), 'Color', [228,26,28]/255, 'LineWidth', 2)
errorbar(1:Nsess, lapsemean(1:Nsess), lapseerr(1:Nsess), 'o', ...,
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
hold on
plot(xvals, ypred, 'LineWidth', 2, 'Color', bluecol);

mymakeaxis('x_label', 'Session', 'y_label', 'Lapse', 'xytitle', 'N = 15')

%%

function loss = lossfun(p, x, y)
A = p(1);
k = p(2);
c = p(3);
ypred = A * exp(-k*x) + c;
loss = sum((ypred - y).^2);
end


function loss = losssigmoid(p, x, y)
lapseL = 0.5;
slope = p(1);
offset = p(2);
lapseR = p(3);
ypred = lapseL + (1 - lapseL - lapseR) ./ (1 + exp(-slope * (x - offset)));
loss = sum((ypred - y).^2);

end