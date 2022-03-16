%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')
paths = pathsetup('matchingsim');
version = '121021';


for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(efflist, prewlst, pswitchlst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
        'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
    plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
    plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

    savefilename = fullfile(paths.figpath, 'infbasedFigs', version,...
        sprintf('inf-based-eff_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end



%% Plotting the switch offset
for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(-PLoffsetlist, prewlst, pswitchlst, 'clim', [0, 10], ...
    'legendname', 'Offset', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
    plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
    plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

    savefilename = fullfile(paths.figpath, 'infbasedFigs', version,...
        sprintf('inf-based-offset_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end


%% Plotting the switch slope


for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(PLslopelist, prewlst, pswitchlst, 'clim', [0, 15], 'legendname', 'Slope', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
    plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
    plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

    savefilename = fullfile(paths.figpath, 'infbasedFigs', version,...
        sprintf('inf-based-slope_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end

%% Plotting the lapse rate (exploration)

for prob = 0:0.1:0.4
    filedir = fullfile(paths.simdatapath, ...
        version, sprintf('EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', ...
        prob, 1 - prob, version));
    load(filedir);
    assert(nblocks == 1000);

    produce_heatmap(LapseR, prewlst, pswitchlst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
    plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
    plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

    savefilename = fullfile(paths.figpath, 'infbasedFigs', version,...
        sprintf('inf-based-lapse_%s_prob%.2f.pdf', version, prob));
    ioutils.savesafe(savefilename, 'pdf'); 
end


