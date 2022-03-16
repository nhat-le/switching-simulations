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
%     ioutils.savesafe(savefilename, 'pdf'); 
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