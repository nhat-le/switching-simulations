rng(124);

plottype = 1; % 1 for constant eps
paperaesthetics;

% For simulating constant eps, vary gamma, eps = 0.2
if plottype
    [offset1, slope1, lapse1] = get_behavioral_params(1, 1);
    [offset2, slope2, lapse2] = get_behavioral_params(3, 3);
    [offset3, slope3, lapse3] = get_behavioral_params(6, 6);
else
% For simulating constant gamma, vary eps, gamma = 1.2
    [offset1, slope1, lapse1] = get_behavioral_params(22, 1);
    [offset2, slope2, lapse2] = get_behavioral_params(22, 9);
    [offset3, slope3, lapse3] = get_behavioral_params(22, 16);
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


% mymakeaxis('x_label', 'Trials from switch', ...
%     'yticks', [1 2 3], 'yticklabels', {'\gamma = 0.1', '\gamma = 0.5', '\gamma = 1.2'})


mymakeaxis('x_label', 'Trials from switch', ...
    'yticks', [1 2 3], 'yticklabels', {'P_s = 0.01, P_r = 0.55', 'P_s = 0.2, P_r = 0.7', 'P_s = 0.45, P_r = 0.99'})





function [offset, slope, lapse] = get_behavioral_params(gammaid, epsid)
% Behavioral params
filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-prob0.0to1.0-071621.mat';
load(filedir);
offset = PRoffsetlist(gammaid, epsid);
slope = PRslopelist(gammaid, epsid);
lapse = LapseR(gammaid, epsid);
end
