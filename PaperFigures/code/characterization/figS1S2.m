%% Plotting the performance
addpath('/Users/minhnhatle/Dropbox (MIT)/Jazayeri/NoisyMutualInhibition/PlotTools')
prob = 0;
type = 'infbased';

expdate = '121021';
rootdir = fullfile('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/', expdate);

if strcmp(type, 'qlearning')
    savedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/qlearningFigs';
%     filedir = sprintf('%sEGreedyQLearningAgent-withCorr-prob%.2fto%.2f-072321.mat', rootdir, prob, 1-prob);  
    filedir = sprintf('%s/EGreedyQLearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', rootdir, prob, 1-prob, expdate);  
elseif strcmp(type, 'infbased')
    savedir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/infbasedFigs';
%     filedir = sprintf('%sEGreedyInferenceBasedAgent-withCorr-prob%.2fto%.2f-072121.mat', rootdir, prob, 1-prob);
    filedir = sprintf('%s/EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat', rootdir, prob, 1-prob, expdate);
end
savedir = sprintf('%s/prob%.2f', savedir, prob);
if ~exist(savedir, 'dir')
    mkdir(savedir);
end


load(filedir);


if strcmp(type, 'qlearning')
produce_heatmap(efflist, epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);
elseif strcmp(type, 'inf-based')
    
produce_heatmap(efflist, prewlst, pswitchlst, 'clim', [0.5,1], 'legendname', 'Efficiency', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

end


% saveas(gcf, sprintf('%s/%s-effprob%.2f.pdf', savedir, type, prob))


%% Plotting the switch offset
if strcmp(type, 'qlearning')
produce_heatmap(-PLoffsetlist, epslst, gammalst, 'clim', [0 10], 'legendname', 'Offset', ...
    'x_label', '$\epsilon$','y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);
elseif strcmp(type, 'inf-based')
produce_heatmap(-PLoffsetlist, prewlst, pswitchlst, 'clim', [0, 10], ...
    'legendname', 'Offset', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);    
    
end
% saveas(gcf, sprintf('%s/%s-offsetprob%.2f.pdf', savedir, type, prob))


%% Plotting the switch slope
if strcmp(type, 'qlearning')
produce_heatmap(PLslopelist, epslst, gammalst, 'clim', [0, 5], 'legendname', 'Slope', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);
elseif strcmp(type, 'inf-based')
produce_heatmap(PLslopelist, prewlst, pswitchlst, 'clim', [0, 15], 'legendname', 'Slope', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);  
end
% saveas(gcf, sprintf('%s/%s-effslope%.2f.pdf', savedir, type, prob))


%% Plotting the lapse rate (exploration)
if strcmp(type, 'qlearning')
produce_heatmap(LapseR, epslst, gammalst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);
elseif strcmp(type, 'inf-based')
produce_heatmap(LapseR, prewlst, pswitchlst, 'clim', [0, 0.5], 'legendname', 'Lapse', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);
plot([0.55, 0.99], [0.01, 0.45], 'k--', 'LineWidth', 2)
plot([0.55, 0.7456, 0.99], [0.01, 0.1986, 0.45], 'kx', 'MarkerSize', 10, 'LineWidth', 2);

  
end
% saveas(gcf, sprintf('%s/%s-efflapse%.2f.pdf', savedir, type, prob))


%% Correlations
% produce_heatmap(ParamsA, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrA', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');
% produce_heatmap(ParamsB, epslst, gammalst, 'clim', [0,10], 'legendname', 'CorrB', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');

slope_at_X4 = ParamsA .* ParamsB .* exp(-ParamsA * 4);

if strcmp(type, 'qlearning')
produce_heatmap(slope_at_X4, epslst, gammalst, 'clim', [0,3], 'legendname', 'CorrB', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'vertline', 0.24, 'horline', 1.22);
elseif strcmp(type, 'infbased')
slope_at_X4 = ParamsA .* ParamsB .* exp(-ParamsA * 4);
produce_heatmap(slope_at_X4, prewlst, pswitchlst, 'clim', [0,3], 'legendname', 'CorrA', ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4); 
end
% saveas(gcf, sprintf('%s/%s-effCorrSlopeX4%.2f.pdf', savedir, type, prob))

% close all
% clear

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


% close all
% clear







