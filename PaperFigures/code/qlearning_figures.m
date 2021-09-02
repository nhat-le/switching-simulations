%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')


filedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat';
load(filedir);

produce_heatmap(efflist, epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);
% plot([epslst(10) epslst(10) epslst(10)], [gammalst(1) gammalst(10) gammalst(22)],  'x', ...
%     'MarkerSize', 9, 'LineWidth', 1)
% plot([epslst(1) epslst(9) epslst(16)], [gammalst(22) gammalst(22) gammalst(22)],  'ko', ...
%     'MarkerSize', 9, 'LineWidth', 1)

%%
meanslope = mean(Qslope_arr, 3);
produce_heatmap(meanslope, epslst, gammalst, 'clim', [0 5], 'legendname', 'Efficiency', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);

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
    'x_label', '$\epsilon$','y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);


%% Plotting the switch slope
produce_heatmap(PLslopelist, epslst, gammalst, 'clim', [0, 5], 'legendname', 'Slope', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);

%% Plotting the lapse rate (exploration)
produce_heatmap(LapseR, epslst, gammalst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);


%% Correlations
% produce_heatmap(ParamsA, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrA', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');
% produce_heatmap(ParamsB, epslst, gammalst, 'clim', [0,10], 'legendname', 'CorrB', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');

slope_at_X4 = ParamsA .* ParamsB .* exp(-ParamsA * 4);
produce_heatmap(slope_at_X4, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrB', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);


%% Theoretical eff
% theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
% produce_heatmap(theoryeff, epslst, gammalst, 'clim', [0.5 1], 'Theory eff');


%% Plot correlation fits
pAarr = ParamsA(:,12);
pBarr = ParamsB(:,12);
xarr = 1:20;
colors = linspace(0.3, 1, 5);

figure;
lines = [];
labels = {};
for id = 1:5
%     subplot(2,4,id)
    i = 1 + (id - 1) * 4;
    A = pAarr(i);
    B = pBarr(i);
    yarr = B - B*exp(-A*xarr);
    lines(id) = plot(xarr, yarr, 'Color', [0 0 colors(id)]);
    labels{id} = sprintf('%.2f', gammalst(i));
    hold on
    ylim([0, 10])
%     title(sprintf('gamma = %.2f', gammalst(i)));
end

mymakeaxis('x_label', 'Number of rewards', 'y_label', 'Number of errors');
c = legend(lines, labels);
c.Title.String = '\gamma';







