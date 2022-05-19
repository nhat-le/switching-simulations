%% rho variations on Q-learning and inference-based spaces (Fig. 2d)

load('data/corr_arr_simulations.mat')

kernel = ones(8) / 64;
imsmoothQ = conv2(corr_arr_Q, kernel, 'same');
imsmoothIB = conv2(corr_arr_ib, kernel, 'same');


produce_heatmap(corr_arr_Q, eps, gamma, 'clim', [0,0.5], 'legendname', '\rho', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', 'ytickvalues', 0:0.2:0.8,...
    'font_size', 22);
hold on
plot([0.1], [0.6], 'k*', 'LineWidth', 1, 'MarkerSize', 30)
plot([0.01], [0.1], 'kx', 'LineWidth', 1, 'MarkerSize', 30)


produce_heatmap(corr_arr_ib(1:end-1,:), prew, psw, 'clim', [0,0.5], 'legendname', '\rho', ...
    'x_label', '$P_{reward}$', 'y_label', '$P_{switch}$', 'ytickvalues', 0:0.1:0.5,...
    'font_size', 22);
hold on
plot([0.7], [0.1],  'kx', 'LineWidth', 1, 'MarkerSize', 30)
plot([0.55], [0.02],  'k*', 'LineWidth', 1, 'MarkerSize', 30)



%% Plot the correlations (Fig. 2a,b right panels + Fig. 2c)
load('data/infbased_agent1_rew_err_correlations.mat')

% scatter(Nr, Ne, 'o', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', 'b')
figure;
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


%% Plot the exemplar agent Q-learning (Fig. 2a left panel)
load data/rew_err_correlations_schematic_qlearning.mat

% figure('Position', [440,487,679,311]);
targets = mod(0:size(choicelst, 1)-1, 2);
targetlst = repmat(targets', [1 size(choicelst, 2)]);

blockrange = 3:15;

q1narrow = q1lst(blockrange,:);
q1flat = reshape(q1narrow', [], 1);
q1flat = q1flat(~isnan(q1flat));

q0narrow = q0lst(blockrange,:);
q0flat = reshape(q0narrow', [], 1);
q0flat = q0flat(~isnan(q0flat));

subtargets = targetlst(blockrange, :);
subtargets = reshape(subtargets', [], 1);

outcomes_narrow = outcomelst(blockrange,:);
outcomes_flat = reshape(outcomes_narrow', [], 1);
outcomes_flat = outcomes_flat(~isnan(outcomes_flat));

choices_narrow = choicelst(blockrange,:);
choices_flat = reshape(choices_narrow', [], 1);
subtargets = subtargets(~isnan(choices_flat));
choices_flat = choices_flat(~isnan(choices_flat));

outcomes = outcomes_flat;
idx = 1:numel(outcomes);

cols = paperaesthetics;

figure('Position', [372,277,936,393]);
l1 = plot(q1flat, 'b', 'LineWidth', 0.5);
hold on
l2 = plot(q0flat, 'k', 'LineWidth', 0.5);

plot(idx(outcomes == 1), choices_flat(outcomes == 1), 'o', ...
    'MarkerFaceColor', cols.bluecol, 'MarkerEdgeColor', 'w', 'MarkerSize', 10)
plot(idx(outcomes == 0), choices_flat(outcomes == 0), 'x', ...
    'Color', cols.redcol, 'MarkerSize', 10)


%plot vertical transition lines
transitions = find(diff(subtargets));
vline(transitions + 0.5);
xlim([0 170])

mymakeaxis('x_label', 'Trials', 'y_label', 'Choice or Value', 'xticks', 0:50:150)
legend([l1, l2], {'q_R', 'q_L'}, 'FontSize', 16);

%% Plot the exemplar agent inf-based (Fig. 2b left panel)
load data/rew_err_correlations_schematic_infbased.mat

% figure('Position', [440,487,679,311]);
targets = mod(0:size(choicelst, 1)-1, 2);
targetlst = repmat(targets', [1 size(choicelst, 2)]);

blockrange = 10:22;

subtargets = targetlst(blockrange, :);
subtargets = reshape(subtargets', [], 1);

choices_narrow = choicelst(blockrange,:);
choices_flat = reshape(choices_narrow', [], 1);

outcomes_narrow = outcomelst(blockrange,:);
outcomes_flat = reshape(outcomes_narrow', [], 1);
outcomes_flat = outcomes_flat(~isnan(outcomes_flat));


subtargets = subtargets(~isnan(choices_flat));
choices_flat = choices_flat(~isnan(choices_flat));

outcomes = outcomes_flat; %subtargets == choices_flat;
idx = 1:numel(outcomes);

cols = paperaesthetics;

figure('Position', [372,277,936,393]);
% l1 = plot(q1flat, 'b', 'LineWidth', 0.5);
hold on
% l2 = plot(q0flat, 'k', 'LineWidth', 0.5);

plot(idx(outcomes == 1), choices_flat(outcomes == 1), 'o', ...
    'MarkerFaceColor', cols.bluecol, 'MarkerEdgeColor', 'w', 'MarkerSize', 10)
plot(idx(outcomes == 0), choices_flat(outcomes == 0), 'x', ...
    'Color', cols.redcol, 'MarkerSize', 10)


%plot vertical transition lines
transitions = find(diff(subtargets));
vline(transitions + 0.5);
xlim([0 150])

mymakeaxis('x_label', 'Trials', 'y_label', 'Choice or Value', 'xticks', 0:50:250)






