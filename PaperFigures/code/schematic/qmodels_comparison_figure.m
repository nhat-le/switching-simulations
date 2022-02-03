paths = pathsetup('matchingsim');

opts.filedir = fullfile(paths.simdatapath, '07232021/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat');
load(opts.filedir);
rng(124);

paperaesthetics;

plottype = 0; % 1 for constant eps

% For simulating constant eps, vary gamma, eps = 0.2
if plottype
    [offset1, slope1, lapse1] = get_behavioral_params(1, 10, opts);
    [offset2, slope2, lapse2] = get_behavioral_params(10, 10, opts);
    [offset3, slope3, lapse3] = get_behavioral_params(22, 10, opts);
else
% For simulating constant gamma, vary eps, gamma = 1.2
    [offset1, slope1, lapse1] = get_behavioral_params(22, 1, opts);
    [offset2, slope2, lapse2] = get_behavioral_params(22, 9, opts);
    [offset3, slope3, lapse3] = get_behavioral_params(22, 16, opts);
end

Nsteps = 25;

% generate behavior
xvals = 1:Nsteps;
prob1 = lapse1 + (1 - 2*lapse1) ./ (1 + exp(slope1 * (xvals + offset1)));
prob2 = lapse2 + (1 - 2*lapse2) ./ (1 + exp(slope2 * (xvals + offset2)));
prob3 = lapse2 + (1 - 2*lapse3) ./ (1 + exp(slope3 * (xvals + offset3)));

if ~plottype
    choice1 = zeros(1, numel(xvals) + 1);
else
    choice1 = [0 rand(1, numel(xvals)) < prob1];
end
choice2 = [0 rand(1, numel(xvals)) < prob2];
choice3 = [0 rand(1, numel(xvals)) < prob3];


%% Show
cols = paperaesthetics;
redcol = cols.redcol;
bluecol = cols.bluecol;

xvals = 0:Nsteps;
if plottype
    plot(xvals(choice1 == 1), 1, 'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol,...
        'MarkerSize', 10);
end
hold on
plot(xvals(choice1 == 0), 1, 'x', 'MarkerSize', 14, 'LineWidth', 2, 'Color', redcol);
plot(xvals(choice2 == 1), 2, 'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol,...
    'MarkerSize', 10);
plot(xvals(choice2 == 0), 2, 'x', 'MarkerSize', 14, 'LineWidth', 2, 'Color', redcol);
plot(xvals(choice3 == 1), 3, 'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol,...
    'MarkerSize', 10);
plot(xvals(choice3 == 0), 3, 'x', 'MarkerSize', 14, 'LineWidth', 2, 'Color', redcol);

vline(0, 'k--')


if plottype
    mymakeaxis('x_label', 'Trials from switch', ...
        'yticks', [1 2 3], 'yticklabels', {'\gamma = 0.1', '\gamma = 0.5', '\gamma = 1.2'})
else
    mymakeaxis('x_label', 'Trials from switch', ...
        'yticks', [1 2 3], 'yticklabels', {'\epsilon = 0.01', '\epsilon = 0.2', '\epsilon = 0.5'})
end




function [offset, slope, lapse] = get_behavioral_params(gammaid, epsid, opts)
% Behavioral params
filedir = opts.filedir;
load(filedir);
offset = PRoffsetlist(gammaid, epsid);
slope = PRslopelist(gammaid, epsid);
lapse = LapseR(gammaid, epsid);
end
