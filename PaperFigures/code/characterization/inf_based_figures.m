%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')

filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/12102021/EGreedyinf-basedAgent-withCorr-doublesigmoid-prob0.30to0.70-121021.mat';
% filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-withCorr-prob0.00to1.00-072121.mat';
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
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

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
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);


%% Plotting the switch slope
produce_heatmap(PLslopelist, prewlst, pswitchlst, 'clim', [0, 15], 'legendname', 'Slope', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

%% Plotting the lapse rate (exploration)
produce_heatmap(LapseR, prewlst, pswitchlst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);


%% Correlations
% produce_heatmap(ParamsA, prewlst, pswitchlst, 'clim', [0,3], 'legendname', 'CorrA', ...
%     'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
% produce_heatmap(ParamsB, prewlst, pswitchlst, 'clim', [0,10], 'legendname', 'CorrB', ...
%     'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
slope_at_X4 = ParamsA .* ParamsB .* exp(-ParamsA * 4);
produce_heatmap(slope_at_X4, prewlst, pswitchlst, 'clim', [0,3], 'legendname', 'CorrA', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);

%% Theoretical eff
% theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
% produce_heatmap(theoryeff, prewlst, pswitchlst, [0.5 1], 'Theory eff');

% a(1:5)
% 
% a(0:5)
% 
% 
% 0 , 5, 10 , 15
% i* 5

%% Plot correlation fits
pAarr = ParamsA(:,2);
pBarr = ParamsB(:,2);
xarr = 1:20;
colors = linspace(0.3, 1, 5);

figure;
lines = [];
labels = {};
plotlst = [1,2,3,5, 10];
for id = 1:5
%     subplot(2,5,id)
    i = plotlst(id);
    A = pAarr(i);
    B = pBarr(i);
    yarr = B - B*exp(-A*xarr);
    lines(id) = plot(xarr, yarr, 'Color', [0 0 colors(id)]);
    labels{id} = sprintf('%.2f', prewlst(i));
    hold on
    ylim([0, 10])
%     title(sprintf('gamma = %.2f', gammalst(i)));
end

mymakeaxis('x_label', 'Number of rewards', 'y_label', 'Number of errors');
c = legend(lines, labels);
c.Title.String = 'P_r';
