%% Diagram for q-learning agent
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/qlearning-sim2.mat')
figure('Position', [0 0 450 250]);
trialid = 1:numel(choices);
plot(trialid(rewards == 1 & choices == 1), choices(rewards == 1 & choices == 1), 'ro', 'MarkerFaceColor', 'r')

hold on
plot(trialid(rewards == 1 & choices == 0), choices(rewards == 1 & choices == 0), 'bo', 'MarkerFaceColor', 'b')
plot(trialid(rewards == 0 & choices == 1), choices(rewards == 0 & choices == 1), 'rx', 'MarkerFaceColor', 'r')
plot(trialid(rewards == 0 & choices == 0), choices(rewards == 0 & choices == 0), 'bx', 'MarkerFaceColor', 'b')


l = vline(31, 'k--');
set(l, 'LineWidth', 1.5);
l1 = plot(q0, 'b', 'LineWidth', 1.5);
l2 = plot(q1, 'k', 'LineWidth', 1.5);
legend([l1 l2], {'q_0', 'q_1'})
title('Q-learning, \gamma = 0.4')

xlabel('Trials')
ylabel('Choice or Probability')
set(gca, 'FontSize', 16)
set(gca, 'box', 'off')
set(gca, 'color', 'none')


%% Matching agent
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/matching-sim.mat');
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
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/matching-sim2.mat')
choices(1) = 0;
choices(2:32) = 1;
choices(33:end) = 0;

rewards(1) = 0;
rewards(2:32) = 1;
rewards(32) = 0;
rewards(33:end) = 1;

prob(1) = 0.5;
prob(2:32) = 0.99;
prob(33:60) = 0.01;



figure('Position', [0 0 450 250]);
trialid = 1:numel(choices);
plot(trialid(rewards == 1 & choices == 1), choices(rewards == 1 & choices == 1), 'ro', 'MarkerFaceColor', 'r')

hold on
plot(trialid(rewards == 1 & choices == 0), choices(rewards == 1 & choices == 0), 'bo', 'MarkerFaceColor', 'b')
plot(trialid(rewards == 0 & choices == 1), choices(rewards == 0 & choices == 1), 'rx', 'MarkerFaceColor', 'r')
plot(trialid(rewards == 0 & choices == 0), choices(rewards == 0 & choices == 0), 'bx', 'MarkerFaceColor', 'b')


prob = prob * 0.8 + 0.1;

l = vline(31, 'k--');
set(l, 'LineWidth', 1.5);
l2 = plot(prob, 'k', 'LineWidth', 1.5);
title('Model-based inference')

xlabel('Trials')
ylabel('Choice or Probability')
set(gca, 'FontSize', 16)
set(gca, 'box', 'off')
set(gca, 'color', 'none')


% targets = [0 0 0 0 0 0 1 1 1 1 1 1];
