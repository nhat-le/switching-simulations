%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')


filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat';
load(filedir);

produce_heatmap(efflist, epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');

%%
% eff_smallgamma = PLslopelist(gammalst < 1, :);
% [gammamax,epsmax] = find(eff_smallgamma == max(eff_smallgamma(:)));
% gamma_maxval = gammalst(gammamax);
% eps_maxval = epslst(epsmax);
% plot(eps_maxval, gamma_maxval, 'w+', 'MarkerSize', 20)
% imshow(efflist2, 'InitialMagnification', 'fit', 'XData', [min(xlst), max(xlst)],...
%     'YData', [min(ylst), max(ylst)]);
% colormap hot
% set(gca, 'YDir', 'normal');
% set(gca, 'XDir', 'normal');
% set(gca, 'Visible', 'on')
% axis square

% What is the optimal value for gamma < 1?


%% Plotting the switch offset
produce_heatmap(-PLoffsetlist, epslst, gammalst, 'clim', [0 10], 'legendname', 'Offset', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');


%% Plotting the switch slope
produce_heatmap(PLslopelist, epslst, gammalst, 'clim', [0, 5], 'legendname', 'Slope', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$')

%% Plotting the lapse rate (exploration)
produce_heatmap(LapseR, epslst, gammalst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');


%% Correlations
produce_heatmap(ParamsA, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrA', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');
produce_heatmap(ParamsB, epslst, gammalst, 'clim', [0,10], 'legendname', 'CorrB', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$');


%% Theoretical eff
% theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
% produce_heatmap(theoryeff, epslst, gammalst, 'clim', [0.5 1], 'Theory eff');
