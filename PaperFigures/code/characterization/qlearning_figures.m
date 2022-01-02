%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')

paths = pathsetup('matchingsim');
version = '121021';


for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyqlearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(efflist, epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.1, 'horline', 1.2);

    savefilename = fullfile(paths.figpath, 'qlearningFigs', version,...
        sprintf('qlearning-eff_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end


%% Plotting the switch offset

for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyqlearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(-PLoffsetlist, epslst, gammalst, 'clim', [0 10], 'legendname', 'Offset', ...
    'x_label', '$\epsilon$','y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);

    savefilename = fullfile(paths.figpath, 'qlearningFigs', version,...
        sprintf('qlearning-offset_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end


%% Plotting the switch slope

for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyqlearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(PLslopelist, epslst, gammalst, 'clim', [0, 5], 'legendname', 'Slope', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);

    savefilename = fullfile(paths.figpath, 'qlearningFigs', version,...
        sprintf('qlearning-slope_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end


%% Plotting the lapse rate (exploration)

for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyqlearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(LapseR, epslst, gammalst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);

    savefilename = fullfile(paths.figpath, 'qlearningFigs', version,...
        sprintf('qlearning-lapse_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end


%% Deprecated code
% Correlations
% produce_heatmap(ParamsA, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrA', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');
% produce_heatmap(ParamsB, epslst, gammalst, 'clim', [0,10], 'legendname', 'CorrB', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');

% slope_at_X4 = ParamsA .* ParamsB .* exp(-ParamsA * 4);
% produce_heatmap(slope_at_X4, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrB', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);


%% Theoretical eff
% theoryeff = 0.5 + (1 - LapseL) .* (0.5 + PLoffsetlist / 30);
% produce_heatmap(theoryeff, epslst, gammalst, 'clim', [0.5 1], 'Theory eff');


%% Plot correlation fits
% pAarr = ParamsA(:,12);
% pBarr = ParamsB(:,12);
% xarr = 1:20;
% colors = linspace(0.3, 1, 5);
% 
% figure;
% lines = [];
% labels = {};
% for id = 1:5
% %     subplot(2,4,id)
%     i = 1 + (id - 1) * 4;
%     A = pAarr(i);
%     B = pBarr(i);
%     yarr = B - B*exp(-A*xarr);
%     lines(id) = plot(xarr, yarr, 'Color', [0 0 colors(id)]);
%     labels{id} = sprintf('%.2f', gammalst(i));
%     hold on
%     ylim([0, 10])
% %     title(sprintf('gamma = %.2f', gammalst(i)));
% end
% 
% mymakeaxis('x_label', 'Number of rewards', 'y_label', 'Number of errors');
% c = legend(lines, labels);
% c.Title.String = '\gamma';



%%
% meanslope = mean(Qslope_arr, 3);
% produce_heatmap(meanslope, epslst, gammalst, 'clim', [0 5], 'legendname', 'Efficiency', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);

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




