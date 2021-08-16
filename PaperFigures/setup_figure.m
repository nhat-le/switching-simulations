%% Figure of the features of switching
% Generate behavior according to a sigmoidal with lapse
xlst = -10:40;

rng(125)
lapse = 0.1;
offset = -12;
slope = 1;
nblockstart = 0; % Block starts on this trial
nblockend = 30;
trialside = (xlst >= nblockstart) & (xlst < nblockend); % Simulating a block switch at 0

rates = lapse + (1 - 2 * lapse) * 1 ./ (1 + exp(-slope * (xlst + offset)));

% Generate behavior based on rates
choices = rand(size(rates)) < rates;

corrtrials = choices == trialside;

% Plot
figure('Position', [440,317,695,481]);
plot(xlst(corrtrials) - nblockstart, choices(corrtrials), 'o', 'Color', [5,113,176]/255,...
    'MarkerFaceColor', [5,113,176]/255)
hold on
plot(xlst(~corrtrials) - nblockstart, choices(~corrtrials), 'x', 'Color', [202,0,32]/255,...
    'MarkerSize', 9, 'LineWidth', 2)


plot(xlst - nblockstart, rates, 'k', 'LineWidth', 2)
plot([nblockstart nblockstart], [0, 1], 'k--', 'LineWidth', 1.5);
plot([nblockend nblockend], [0, 1], 'k--', 'LineWidth', 1.5);

% Marker for parameters
plot([-offset -offset], [0.2, 0.8], 'k--')
% quiver(20, 0.9, 0, 0.1)
% Lapse annotation
annotation('arrow', 'X', [0.7, 0.7], 'Y', [0.87, 0.92], 'HeadLength', 8, 'HeadWidth', 8);
annotation('arrow', 'X', [0.7, 0.7], 'Y', [0.92, 0.87], 'HeadLength', 8, 'HeadWidth', 8);

% Offset annotation
annotation('arrow', 'X', [0.475, 0.605], 'Y', [0.65, 0.65], 'HeadLength', 8, 'HeadWidth', 8);
annotation('arrow', 'X', [0.605, 0.475], 'Y', [0.65, 0.65], 'HeadLength', 8, 'HeadWidth', 8);


mymakeaxis('x_label', 'Trials in block', 'y_label', 'Choice', 'offsetRatio', 0.05)

%%
text(4.5, 0.45, 'N_S', 'Interpreter', 'tex', 'FontSize', 15, 'FontAngle', 'italic')
text(22, 0.95, 'P_E', 'Interpreter', 'tex', 'FontSize', 15, 'FontAngle', 'italic')
text(14.6, 0.7, '\alpha', 'Interpreter', 'tex', 'FontSize', 15)



