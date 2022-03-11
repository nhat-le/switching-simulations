%% Agent behavior illustration

load shuffle_simulation_results_lapse0.1_013022_v6.mat

% shuffle the sim results
interleave = 1;
Nsubblock = 250;
Nblocks = 1000;
cols = paperaesthetics;
if interleave
    idxall = [];
    % Interleave blocks 1 and 2 in subblocks of 200
    idx = [1:Nsubblock,  (Nblocks + 1) : (Nblocks + Nsubblock)];
    for j = 1:(Nblocks / Nsubblock)
        idxall = [idxall idx];
        idx = idx + Nsubblock;

    end 
    block = block(idxall, :);
end

redbluecolormap = [178,24,43; 209,229,240]/ 255;
statecols = brewermap(4, 'Set1');

figure(1);
subplot(121)
imagesc(1-block)
hold on
colormap(redbluecolormap)
axis xy
ylim([0 2000])

% Plot the regime demarcation
xbar = 16;
for j = 1:(Nblocks * 2 / Nsubblock)
    if mod(j, 2) == 1
        plot([xbar xbar], [(j-1) * Nsubblock + 1 j * Nsubblock + 1], 'Color', statecols(2,:), 'LineWidth', 3)
    else
        plot([xbar xbar], [(j-1) * Nsubblock + 1 j * Nsubblock + 1], 'Color', statecols(1,:), 'LineWidth', 3)       
    end
end

xlim([0 18])
mymakeaxis('x_label', 'Trials in block', 'y_label', 'Block #', 'xytitle', 'Agent 1',...
    'xticks', 0:5:15, 'font_size', 18)



subplot(122)
imagesc(1-shuffleblock)
hold on
colormap(redbluecolormap)
axis xy
ylim([0 2000])
xlim([0 18])
plot([xbar xbar], [1 2*Nblocks], 'Color', statecols(3,:), 'LineWidth', 3)




mymakeaxis('x_label', 'Trials in block', 'y_label', 'Block #', 'xytitle', 'Agent 2',...
    'xticks', 0:5:15, 'font_size', 18)

%% Transition function visualization
load shuffle_simulation_results_lapse0.1_013022_v6.mat
Nblocks = size(block, 1);
statecols = brewermap(4, 'Set1');

figure;
l1 = plot(1-mean(block, 1), 'Color', statecols(4,:), 'LineWidth', 2);
hold on
% plot(1-mean(block(1:Nblocks/2,:), 1))
% plot(1-mean(block(Nblocks/2+1:end,:), 1))
l2 = plot(1-mean(shuffleblock, 1), 'Color', statecols(3,:), 'LineWidth', 2);

ylim([0, 1])

% plot(1 - mean(blockalt, 1))
mymakeaxis('x_label', 'Trials in block', 'y_label', 'P (Correct)', 'font_size', 18,...
    'xticks', 0:5:15)
leg = legend([l1, l2], {'Agent 1', 'Agent 2'}, 'FontSize', 16);


%% Logistic regression results
load shuffle_simulation_results_lapse0.1_013022_v6.mat
N = numel(coef1);
N = N/3;
ymin = -0.4;
ymax = 0.4;
cols = paperaesthetics;
statecols = brewermap(4, 'Set1');

figure('Position', [440,477,843,321]);
subplot(131)
plot(0:N-1, coef1(N:-1:1), 'Color', statecols(4,:), 'LineWidth', 2);
hold on
plot(0:N-1, coef3(N:-1:1), 'Color', statecols(3,:), 'LineWidth', 2);
% set(gca, 'XDir', 'reverse');
ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. choice',...
    'xticks', 0:5:15, 'xticklabels', {'0', '-5', '-10', '-15'}, 'font_size', 18)



subplot(132)
plot(0:N-1, coef1(N*2:-1:N+1), 'Color', statecols(4,:), 'LineWidth', 2)
hold on
plot(0:N-1, coef3(N*2:-1:N+1), 'Color', statecols(3,:), 'LineWidth', 2)
ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. reward',...
    'xticks', 0:5:15, 'xticklabels', {'0', '-5', '-10', '-15'}, 'font_size', 18)


subplot(133)
l1 = plot(0:N-1, coef1(end:-1:N*2+1), 'Color', statecols(4,:), 'LineWidth', 2);
hold on
l2 = plot(0:N-1, coef3(end:-1:N*2+1), 'Color', statecols(3,:), 'LineWidth', 2);
ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. choice x reward',...
    'xticks', 0:5:15, 'xticklabels', {'0', '-5', '-10', '-15'}, 'font_size', 18)
legend([l1, l2], {'Agent 1', 'Agent 2'}, 'FontSize', 18);




%% Q-learning simulation results
load rlmodel_sim_020122_v5.mat
statecols = brewermap(4, 'Set1');
figure('Position', [440,357,808,441])
l = violin([coefslst1', biaslst1'], ...
    'mc','','medc','k' , 'facecolor', statecols(4,:), 'x', [1,3]);
l = violin([coefs_all{1}' bias_all{1}'], ...
    'mc','','medc','k' , 'facecolor', statecols(3,:), 'x', [2,4]);
ylim([-0.05, 0.1])
xlim([0 5])
mymakeaxis('x_label', 'Parameter', 'y_label', ...
    'Coefficient', 'yticks', -0.05:0.05:0.1, 'xticks', [1,2,3,4], ...
    'xticklabels', {'Agent 1', 'Agent 2', 'Agent 1', 'Agent 2'}, 'font_size', 20);
legend('off')


