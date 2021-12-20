% Figures for decoding performance
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/svmresults_from_pickle_090221.mat');
% 
% mindecodingQ = min(decoding_perf, [], [3, 4]);
% mindecodingIB = squeeze(min(decoding_perf, [], [1, 2]));
% 
% %%
% produce_heatmap(mindecodingQ', epslst, gammalst, 'clim', [0.5,1], 'legendname', 'Decoding performance', ...
%     'x_label', '$\epsilon$', 'y_label', '$\gamma$');
% 
% 
% %%
% produce_heatmap(mindecodingIB', prewlst, pswitchlst, 'clim', [0.5, 1], 'legendname', 'Decoding performance', ...
%     'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.4);

