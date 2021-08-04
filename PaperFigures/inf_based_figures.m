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

produce_heatmap_inf_based(efflist, prewlst, pswitchlst, [0.5,1], 'Efficiency', ...
    'xgridsize', 20, 'ygridsize', 25);

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
produce_heatmap_inf_based(-PLoffsetlist, prewlst, pswitchlst, [0 10], 'Offset', 20, 25);


%% Plotting the switch slope
produce_heatmap_inf_based(PLslopelist, prewlst, pswitchlst, [0, 15], 'Slope', 20, 25)

%% Plotting the lapse rate (exploration)
produce_heatmap_inf_based(LapseR, prewlst, pswitchlst, [0, 0.5], 'Lapse', 20, 25);


%% Correlations
produce_heatmap_inf_based(ParamsA, prewlst, pswitchlst, [0,3], 'CorrA', 20, 25);
produce_heatmap_inf_based(ParamsB, prewlst, pswitchlst, [0,10], 'CorrB', 20, 25);


%% Theoretical eff
theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
produce_heatmap(theoryeff, prewlst, pswitchlst, [0.5 1], 'Theory eff');

