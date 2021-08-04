%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')


filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat';
load(filedir);

produce_heatmap(efflist, epslst, gammalst, [0.5,1], 'Efficiency');

%%
eff_smallgamma = PLslopelist(gammalst < 1, :);
[gammamax,epsmax] = find(eff_smallgamma == max(eff_smallgamma(:)));
gamma_maxval = gammalst(gammamax);
eps_maxval = epslst(epsmax);
plot(eps_maxval, gamma_maxval, 'w+', 'MarkerSize', 20)
% imshow(efflist2, 'InitialMagnification', 'fit', 'XData', [min(xlst), max(xlst)],...
%     'YData', [min(ylst), max(ylst)]);
% colormap hot
% set(gca, 'YDir', 'normal');
% set(gca, 'XDir', 'normal');
% set(gca, 'Visible', 'on')
% axis square

% What is the optimal value for gamma < 1?


%% Plotting the switch offset
produce_heatmap(-PLoffsetlist, epslst, gammalst, [0 10], 'Offset');


%% Plotting the switch slope
produce_heatmap(PLslopelist, epslst, gammalst, [0, 5], 'Slope')

%% Plotting the lapse rate (exploration)
produce_heatmap(LapseR, epslst, gammalst, [0, 0.5], 'Lapse');


%% Correlations
produce_heatmap(ParamsA, epslst, gammalst, [0,3], 'CorrA');
produce_heatmap(ParamsB, epslst, gammalst, [0,10], 'CorrB');


%% Theoretical eff
theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
produce_heatmap(theoryeff, epslst, gammalst, [0.5 1], 'Theory eff');
