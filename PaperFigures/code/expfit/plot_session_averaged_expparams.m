% load expdata/expfit_params.mat
paths = pathsetup('matchingsim');
version = '122221b';
filterval = 20; %filter values such that expoffsets_all(expoffsets_all > 20) = filterval;
files = dir(fullfile(paths.expdatapath, version, 'fitparams_session_averaged*.mat'));
assert(numel(files) == 1)
load(fullfile(files(1).folder, files(1).name));

% Load the fitrange file
load(fullfile(paths.expdatapath, '102121', 'fitranges_122221.mat'));

Nsess = 30; %plot only this many sessions
cols = paperaesthetics;
bluecol = cols.bluecol;
% bluecol = [8,81,156]/255;


%%
animals_ids = fields(fitparams_all);
effs_all = {};
slopes_all = {};
lapses_all = {};
offsets_all = {};

for i = 1:numel(animals_ids)
    [offsets, slopes, lapses, effs] = parse_params(fitparams_all.(animals_ids{i}));
    
    %Trim based on fitparams
    idx = contains(animals, animals_ids{i});
    assert(sum(idx) == 1);
    fitrange = ranges{idx};
    
    effs_all{i} = effs(fitrange);
    slopes_all{i} = slopes(fitrange);
    lapses_all{i} = lapses(fitrange);
    offsets_all{i} = offsets(fitrange);
end


%% Pad to same length
expeff_all = pad_to_same_length(effs_all);
expoffsets_all = -pad_to_same_length(offsets_all);
expoffsets_all(expoffsets_all > 20) = filterval;


expslopes_all = pad_to_same_length(slopes_all);
explapses_all = pad_to_same_length(lapses_all);



%% Eff
figure;
% plot(expeff_all', 'Color', [55,126,184]/255)
effmeans = nanmean(expeff_all, 1);
efferr = nanstd(expeff_all, [], 1) / sqrt(size(expeff_all, 1));

hold on
% plot(nanmean(expeff_all, 1), 'Color', [228,26,28]/255, 'LineWidth', 2)
errorbar(1:Nsess, effmeans(1:Nsess), efferr(1:Nsess), 'o', ...,
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
mymakeaxis('x_label', 'Session', 'y_label', 'Performance',...
    'yticks', 0.5:0.1:0.7,'font_size', 22)

% Fit the eff mean
xvals = 1:Nsess;
yvals = effmeans(1:Nsess);
p = [1, 10, 0.2, 0.5];
pfit = fminsearch(@(p) losssigmoid(p, xvals, yvals), p);
ypred = pfit(4) + (1 - pfit(4) - pfit(3)) ./ (1 + exp(-pfit(1) * (xvals - pfit(2))));
ypred2 = pfit(4) + (1 - pfit(4) - p(3)) ./ (1 + exp(-p(1) * (xvals - p(2))));

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
mymakeaxis('x_label', 'Session', 'y_label', 'Offset',...
    'font_size', 22)


%% Slopes
figure;
% plot(expeff_all', 'Color', [55,126,184]/255)
slopemeans = nanmean(expslopes_all, 1);
slopeerr = nanstd(expslopes_all, [], 1) / sqrt(size(expslopes_all, 1));

hold on
% plot(nanmean(expeff_all, 1), 'Color', [228,26,28]/255, 'LineWidth', 2)
errorbar(1:Nsess, slopemeans(1:Nsess), slopeerr(1:Nsess), 'o', ...,
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
mymakeaxis('x_label', 'Session', 'y_label', 'Slope',...
    'font_size', 22)

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

mymakeaxis('x_label', 'Session', 'y_label', 'Lapse','font_size', 22)

%% Save the exp metrics
filename = fullfile(paths.expdatapath, version, sprintf('expfit_params_%s.mat', version));
if ~exist(filename, 'file')
    save(filename, 'expeff_all', 'explapses_all', 'expslopes_all', 'expoffsets_all',...
        'animals_ids', 'Nsess', 'filterval');
    fprintf('File saved!\n')
else
    fprintf('File exists, skipping save...\n');
end


function [offsets, slopes, lapses, effs] = parse_params(params)

slopesL = params.pL(:,1);
slopesR = params.pR(:,1);
offsetsL = params.pL(:,2);
offsetsR = params.pR(:,2);
lapsesL = params.pL(:,3);
lapsesR = params.pR(:,3);

slopes = nanmin([slopesL, slopesR], [], 2)';
offsets = nanmean([offsetsL, offsetsR], 2)';
lapses = nanmax([lapsesL, lapsesR], [], 2)';
effs  = params.eff;




end


function loss = lossfun(p, x, y)
A = p(1);
k = p(2);
c = p(3);
ypred = A * exp(-k*x) + c;
loss = sum((ypred - y).^2);
end


function loss = losssigmoid(p, x, y)
lapseL = p(4);
slope = p(1);
offset = p(2);
lapseR = p(3);
ypred = lapseL + (1 - lapseL - lapseR) ./ (1 + exp(-slope * (x - offset)));
loss = sum((ypred - y).^2);

end



