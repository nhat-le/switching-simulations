%% Q-learning space
path = pathsetup('matchingsim');
paperaesthetics;
load(fullfile(path.simdatapath, 'deprecated/092321', ...
    'EGreedyinf-basedAgent-withCorr-doublesigmoid-prob0.00to1.00-092321.mat'),...
    'epslst', 'gammalst', 'prewlst', 'pswitchlst');

%Grids
[eps, gamma] = meshgrid(epslst, gammalst);
[prew, psw] = meshgrid(prewlst, pswitchlst);

%Q-learning space
cols = paperaesthetics;
bluecol = cols.bluecol;
redcol = cols.redcol;
figure('Position', [440,116,427,682]);
subplot(211)
plot(eps, gamma, 'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol,...
    'MarkerSize', 2)
mymakeaxis('x_label', '\epsilon', 'y_label', '\gamma', 'xytitle', 'Q-learning')



subplot(212)
plot(prew, psw, 'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol,...
    'MarkerSize', 2)
mymakeaxis('x_label', 'P_{rew}', 'y_label', 'P_{switch}', 'xytitle', 'Inference-based')


%% Parameter space colored
% opts = struct;
% [res1, opts] = load_and_run(0, opts);











