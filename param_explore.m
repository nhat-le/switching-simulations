%% Q-learning
load EGreedyQLearning-prob0.2to0.8-071121.mat

PLoffsetQ = PLoffsetlist;
PLslopeQ = PLslopelist;
LapseQ = LapseL;
figure;
imagesc(LapseQ, 'XData', [min(epslst) max(epslst)], 'YData', [min(gammalst), max(gammalst)])
axis xy
% xticks(1:20)
% yticks(1:30)
% xticklabels(round(epslst, 1))
% yticklabels(round(gammalst, 1))

% set(gca, 'XData', [-10 10]);


%% Inference-based
load EGreedyInferenceBasedAgent-prob0.2to0.8-071121.mat
PLoffsetIB = PLoffsetlist;
PLslopeIB = PLslopelist;
LapseIB = LapseL;
figure;
imagesc(LapseL, 'XData', [min(prewlst) max(prewlst)], 'YData', [min(pswitchlst), max(pswitchlst)])
axis xy


%% Compare the two
figure;
subplot(131)
plot(PLoffsetQ(10,:));
hold on
plot(PLoffsetIB(10,:));

subplot(132)
plot(PLslopeQ(10,:));
hold on
plot(PLslopeIB(10,:));


subplot(133)
plot(LapseQ(10,:));
hold on
plot(LapseIB(10,:));


%%
figure;
plot3(PLoffsetIB(:), PLslopeIB(:), LapseIB(:), 'bo')
hold on
plot3(PLoffsetQ(:), PLslopeQ(:), LapseQ(:), 'x')


%%
offsetIBflat = reshape(PLoffsetIB, 1, []);
slopeIBflat = reshape(PLslopeIB, 1, []);
[xx, yy] = meshgrid(prewlst, pswitchlst);
xxflat = reshape(xx, 1, []);
yyflat = reshape(yy, 1, []);

% locs = (offsetIBflat < -1.65) & (slopeIBflat > -22) & (slopeIBflat < -14);
% locs = (offsetIBflat > -1.3) & (slopeIBflat > -22) & (slopeIBflat < -14);
locs = slopeIBflat > -10;
xxgood = xxflat(locs);
yygood = yyflat(locs);
plot(xxgood, yygood, 'o');


%%
offsetQflat = reshape(PLoffsetQ, 1, []);
slopeQflat = reshape(PLslopeQ, 1, []);
[xxQ, yyQ] = meshgrid(epslst, gammalst);
xxflatQ = reshape(xxQ, 1, []);
yyflatQ = reshape(yyQ, 1, []);

% locs = (offsetIBflat < -1.65) & (slopeIBflat > -22) & (slopeIBflat < -14);
% locs = (offsetIBflat > -1.3) & (slopeIBflat > -22) & (slopeIBflat < -14);
locs = (abs(offsetQflat + 1.9) < 0.2) & (abs(slopeQflat + 17) < 4);
xxgood = xxflatQ(locs);
yygood = yyflatQ(locs);
plot(xxgood, yygood, 'o');










