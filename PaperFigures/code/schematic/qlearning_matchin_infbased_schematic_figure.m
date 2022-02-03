%% Diagram for q-learning agent
paths = pathsetup('matchingsim');
paperaesthetics;

load(fullfile(paths.simdatapath, 'schematic/qlearning-sim-122221.mat'))
figure('Position', [0 0 450 250]);
trialid = 1:numel(choices);
plot(trialid(rewards == 1 & choices == 1), choices(rewards == 1 & choices == 1), ...
    'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol, 'Color', bluecol)

hold on
plot(trialid(rewards == 1 & choices == 0), choices(rewards == 1 & choices == 0), 'o', ...
    'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol, 'Color', bluecol)
plot(trialid(rewards == 0 & choices == 1), choices(rewards == 0 & choices == 1), 'x', ...
    'MarkerFaceColor', redcol, 'Color', redcol)
plot(trialid(rewards == 0 & choices == 0), choices(rewards == 0 & choices == 0), 'x', ...
    'MarkerFaceColor', redcol, 'Color', redcol)


l = vline(17, 'k--');
set(l, 'LineWidth', 1.5);
l1 = plot(q0, 'b', 'LineWidth', 1.5);
l2 = plot(q1, 'k', 'LineWidth', 1.5);
legend([l1 l2], {'q_0', 'q_1'})
title(sprintf('Q-learning, \\gamma = %.1f, \\epsilon = %.1f', ...
    params.gammalst, params.epslst))

xlabel('Trials')
ylabel('Choice or Probability')
set(gca, 'FontSize', 16)
set(gca, 'box', 'off')
set(gca, 'color', 'none')


%% Matching agent
load(fullfile(paths.simdatapath, 'schematic/matching-sim.mat'));
figure('Position', [0 0 450 250]);
trialid = 1:numel(choices);
plot(trialid(rewards == 1), choices(rewards == 1), 'o', 'MarkerFaceColor', 'b')

hold on
plot(trialid(rewards == 0), choices(rewards == 0), 'o', 'MarkerFaceColor', 'r')

l = vline(31, 'k--');
set(l, 'LineWidth', 1.5);
l1 = plot(r0/7, 'b', 'LineWidth', 1.5);
l2 = plot(r1/7, 'k', 'LineWidth', 1.5);
legend([l1 l2], {'r_0', 'r_1'})
title('Local matching, \tau = 5')

xlabel('Trials')
ylabel('Choice or Probability')
set(gca, 'FontSize', 16)


%% Inference agent
paths = pathsetup('matchingsim');

load(fullfile(paths.simdatapath, 'schematic/infbased-sim-122221.mat'))
prob = prob1;

figure('Position', [0 0 450 250]);
trialid = 1:numel(choices);
plot(trialid(rewards == 1 & choices == 1), choices(rewards == 1 & choices == 1), ...
    'o', 'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol, 'Color', bluecol)

hold on
plot(trialid(rewards == 1 & choices == 0), choices(rewards == 1 & choices == 0), 'o', ...
    'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol, 'Color', bluecol)
plot(trialid(rewards == 0 & choices == 1), choices(rewards == 0 & choices == 1), 'x', ...
    'MarkerFaceColor', redcol, 'Color', redcol)
plot(trialid(rewards == 0 & choices == 0), choices(rewards == 0 & choices == 0), 'x', ...
    'MarkerFaceColor', redcol, 'Color', redcol)

transition = find(diff(prob0 > 0.5));
l = vline(transition, 'k--');
set(l, 'LineWidth', 1.5);
l2 = plot(prob, 'k', 'LineWidth', 1.5);
title(sprintf('Inference-based agent, p_s = %.1f, p_r = %.1f', ...
    params.pswitchlst, params.prewlst))

xlabel('Trials')
ylabel('Choice or Probability')
set(gca, 'FontSize', 16)
set(gca, 'box', 'off')
set(gca, 'color', 'none')


%% Plot the behavior of the agents over multiple blocks
paths = pathsetup('matchingsim');
load(fullfile(paths.simdatapath, 'schematic/qlearning_sim_combinations.mat'));

figure('Position', [440,368,695,430]);
hold on
colors = brewermap(9, 'Blues');
colors = colors([3, 5, 8], :);
type = 'gamma'; %gamma or epsilon

for id = 1:3
    subplot(2,3,id)
    
    if strcmp(type, 'gamma')
        idchoose = id;
    elseif strcmp(type, 'epsilon')
        idchoose = id + 3;
    else
        error('Invalid idchoose');
    end
    sub_choicearr = squeeze(choicelst_all(idchoose,2:101,:));
    imagesc(1-sub_choicearr, 'AlphaData', ~isnan(sub_choicearr));
    colormap redblue
    axis xy
    mymakeaxis('x_label', 'Trials from switch', 'y_label', 'Block #',...
        'xticks', 0:10:20)
    
    subplot(2,3,id + 3)
    choicelst = squeeze(choicelst_all(idchoose,:,:));

    % imagesc(choicelst)
%     stdshade(choicelst, 0.6, colors(id,:));
    
    % fitted curve
    xvals = -1:15;
    yvals = mathfuncs.sigmoid(xvals, -pfits(idchoose, 1), pfits(idchoose, 2), pfits(idchoose, 3));
    errors = sqrt(nanmean(choicelst) .* (1 - nanmean(choicelst))); % / size(choicelst, 1));
    errorbar(1:size(choicelst, 2), nanmean(choicelst), errors, 'k')
    hold on
    plot(xvals+1, yvals, 'b', 'LineWidth', 2)
    xlim([0, 15])
    ylim([-0.3, 1.3])
    mymakeaxis('x_label', 'Trials from switch', 'y_label', 'P(Correct)', ...
        'xticks', 0:5:15, 'yticks', 0:0.5:1)

    
end



