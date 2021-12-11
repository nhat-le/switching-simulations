%% Reality
% load('/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/f03/2021-02-19/2/2021-02-19_2_F03_Block.mat');
load('/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/f12/2021-06-02/1/2021-06-02_1_F12_Block.mat');
% load('/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/f12/2021-07-08/1/2021-07-08_1_F12_Block.mat');

figure;
ax= axes;
plotPerfSimple(ax, block);
xlim([0, 153])
ylim([-4 4])
mymakeaxis('x_label', 'Trials', 'yticks', [-1, 1])

%% Expectation
load('/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/f03/2021-02-19/2/2021-02-19_2_F03_Block.mat');

block.events.contrastLeftValues = [ones(1, 20) zeros(1, 20) ones(1, 20) zeros(1, 20) ...
    ones(1, 20) zeros(1, 20) 1 1 1] ;
block.events.contrastRightValues =  [zeros(1, 20) ones(1, 20) zeros(1, 20) ones(1, 20) ...
    zeros(1, 20) ones(1, 20) 1 1 1]; 
block.events.responseValues = [1 1 1 ones(1, 20)*-1 ones(1, 20) ones(1, 20)*-1 ones(1, 20)...
    ones(1, 20)*-1 ones(1, 17)];
figure;
ax= axes;
plotPerfSimple(ax, block);
xlim([0, 130])
ylim([-4 4])
mymakeaxis('x_label', 'Trials', 'yticks', [-1, 1])

%% Q learning space
figure;
ax = axes;
xlim([0, 0.5])
ylim([0, 1.5])
mymakeaxis('x_label', '\epsilon', 'y_label', '\gamma')


%% IB space
figure;
ax = axes;
xlim([0.5, 1])
ylim([0, 0.4])
mymakeaxis('x_label', 'P_{rew}', 'y_label', 'P_{switch}')



















