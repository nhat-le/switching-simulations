load corr_arr_simulations_020322.mat

%smooth the Q map
% imsmooth = imfilter(corr_arr_Q);
kernel = ones(8) / 64;
imsmoothQ = conv2(corr_arr_Q, kernel, 'same');
imsmoothIB = conv2(corr_arr_ib, kernel, 'same');


produce_heatmap(corr_arr_Q, eps, gamma, 'clim', [0,0.5], 'legendname', '\rho', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'ytickvalues', 0:0.2:0.8,...
    'font_size', 22);
hold on
plot([0.1], [0.6], 'kx', 'LineWidth', 2, 'MarkerSize', 30)
plot([0.05], [0.05], 'k*', 'LineWidth', 2, 'MarkerSize', 30)


produce_heatmap(corr_arr_ib(1:end-1,:), prew, psw, 'clim', [0,0.5], 'legendname', '\rho', ...
    'x_label', '$P_{reward}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.5,...
    'font_size', 22);
hold on
plot([0.7], [0.1],  'kx', 'LineWidth', 2, 'MarkerSize', 30)
plot([0.55], [0.02],  'k*', 'LineWidth', 2, 'MarkerSize', 30)



%% Plot the exemplar agents
load('Qlearning_agent3_rew_err_correlations.mat')
% load('infbased_agent1_rew_err_correlations.mat')

% scatter(Nr, Ne, 'o', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', 'b')
hold on
errorbar(1:numel(meanvals), meanvals, stdvals, 'o',...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
% ylim([1,3])

cols = paperaesthetics;
% Linear regression
mdl = fitlm(double(Nr), double(Ne));
xvals = 0:max(Nr);
preds = double(xvals) * mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
plot(xvals, preds, 'LineWidth', 2, 'Color', cols.bluecol);
mymakeaxis('x_label', 'Num. rewards', 'y_label', 'Num. errors',...
    'xticks', 0:5:15, 'font_size', 22)

%% Schematic figure
load rew_err_correlations_schematic.mat

figure('Position', [440,487,679,311]);
subplot(121)
q1toplot_b1 = [q1lst(4,:) q1lst(5,:)];
q1toplot_b1 = q1toplot_b1(~isnan(q1toplot_b1));

choicestoplot_b1 = [choicelst(4,:), choicelst(5,:)];
choicestoplot_b1 = choicestoplot_b1(~isnan(choicestoplot_b1));

N1a = sum(~isnan(choicelst(4,:)));
N1b = sum(~isnan(choicelst(5,:)));
targetside = [ones(1, N1a) zeros(1, N1b)];
outcomes = targetside == choicestoplot_b1;
idx = 1:numel(outcomes);

cols = paperaesthetics;
plot(q1toplot_b1)
hold on
plot(idx(outcomes == 1), choicestoplot_b1(outcomes == 1), 'o', ...
    'MarkerFaceColor', cols.bluecol, 'MarkerEdgeColor', cols.bluecol)
plot(idx(outcomes == 0), choicestoplot_b1(outcomes == 0), 'x', ...
    'Color', cols.redcol)

vline(N1a)
xlim([0, 60])
ylim([0, 1])

mymakeaxis('x_label', 'Trials', 'y_label', 'Choice or Value', 'xticks', 0:20:60)

% example block 2
subplot(122)
q1toplot_b2 = [q1lst(8,:) q1lst(9,:)];
q1toplot_b2 = q1toplot_b2(~isnan(q1toplot_b2));

choicestoplot_b2 = [choicelst(8,:), choicelst(9,:)];
choicestoplot_b2 = choicestoplot_b2(~isnan(choicestoplot_b2));

N2a = sum(~isnan(choicelst(8,:)));
N2b = sum(~isnan(choicelst(9,:)));
targetside = [ones(1, N2a) zeros(1, N2b)];
outcomes = targetside == choicestoplot_b2;
idx = 1:numel(outcomes);

cols = paperaesthetics;
plot(q1toplot_b2)
hold on
plot(idx(outcomes == 1), choicestoplot_b2(outcomes == 1), 'o', ...
    'MarkerFaceColor', cols.bluecol, 'MarkerEdgeColor', cols.bluecol)
plot(idx(outcomes == 0), choicestoplot_b2(outcomes == 0), 'x', ...
    'Color', cols.redcol)

vline(N2a)
xlim([0, 60])
ylim([0, 1])

mymakeaxis('x_label', 'Trials', 'y_label', 'Choice or Value', 'xticks', 0:20:60)









