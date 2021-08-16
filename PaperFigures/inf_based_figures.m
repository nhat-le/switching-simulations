%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')


filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/EGreedyInferenceBasedAgent-withCorr-prob0.00to1.00-072121.mat';
load(filedir);

% produce_heatmap(efflist, epslst, gammalst, [0.5,1], 'Efficiency');
%2d interpolation
% [epsgrid, gammagrid] = meshgrid(epslst, gammalst);
% 
% epslstnew = linspace(0.01, 0.5, 80);
% gammalstnew = linspace(0.01, 1.4, 80);
% [epsgridnew, gammagridnew] = meshgrid(epslstnew, gammalstnew);

% xlst = epslst;
% ylst = gammalst;

% efflist2 = interp2(epsgrid,gammagrid,efflist,epsgridnew,gammagridnew);
produce_heatmap(efflist, prewlst, pswitchlst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);

% produce_heatmap(efflist, prewlst, pswitchlst, [0.5,1], 'Efficiency', ...
%     'xgridsize', 20, 'ygridsize', 25);

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
produce_heatmap(-PLoffsetlist, prewlst, pswitchlst, 'clim', [0, 10], ...
    'legendname', 'Offset', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);


%% Plotting the switch slope
produce_heatmap(PLslopelist, prewlst, pswitchlst, 'clim', [0, 15], 'legendname', 'Slope', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);

%% Plotting the lapse rate (exploration)
produce_heatmap(LapseR, prewlst, pswitchlst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);


%% Correlations
produce_heatmap(ParamsA, prewlst, pswitchlst, 'clim', [0,3], 'legendname', 'CorrA', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
produce_heatmap(ParamsB, prewlst, pswitchlst, 'clim', [0,10], 'legendname', 'CorrB', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);


%% Theoretical eff
% theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
% produce_heatmap(theoryeff, prewlst, pswitchlst, [0.5 1], 'Theory eff');

