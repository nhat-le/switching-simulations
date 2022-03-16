%% Data generation
Nblocks = 500;
blocklen = 15;

offset1 = 3;
offset2 = 3;
slope1 = 0.01;
slope2 = 0.8;
lapse1 = 0.1;
lapse2 = 0.1;
Nsubblock = 50;

seed = 124;
rng(seed);
Niter = 10;
Nlogistic = 10;

h = waitbar(0);

gammalst = [];
biaslst = [];
Ball_arr = [];
Ball_arr1 = [];
Ball_arr2 = [];

for i = 1:Niter
    waitbar(i/Niter, h)
    block1 = generate_blocks_with_sigmoid(Nblocks * 2, blocklen, offset1, slope1, lapse1);
    block2 = generate_blocks_with_sigmoid(Nblocks * 2, blocklen, offset2, slope2, lapse2);
    block = [block1(1:Nblocks,:); block2(1:Nblocks,:)];
    idxall = [];
    % Interleave blocks 1 and 2 in subblocks of 200
    idx = [1:Nsubblock,  (Nblocks + 1) : (Nblocks + Nsubblock)];
    for i = 1:(Nblocks / Nsubblock)
        idxall = [idxall idx];
        idx = idx + Nsubblock;

    end 
    block = block(idxall, :);


    %RL fitting
    [choiceslst, outcomeslst] = unfold_block(block);
    res = do_rl_fitting(choiceslst, outcomeslst);

    [choiceslst1, outcomeslst1] = unfold_block(block1);
    res1 = do_rl_fitting(choiceslst1, outcomeslst1);

    [choiceslst2, outcomeslst2] = unfold_block(block2);
    res2 = do_rl_fitting(choiceslst2, outcomeslst2);
    
    gammalst(end+1,:) = [res1(1), res2(1), res(1)];
    biaslst(end+1,:) = [res1(3), res2(3), res(3)];
    
    % Logistic regression
    B_all = do_logistic_regression(block, Nlogistic);
    B_all1 = do_logistic_regression(block1, Nlogistic);
    B_all2 = do_logistic_regression(block2, Nlogistic);
    
    Ball_arr(end+1,:) = B_all;
    Ball_arr1(end+1,:) = B_all1;
    Ball_arr2(end+1,:) = B_all2;

end
close(h)

%% Figure: RL summary
statecols = brewermap(4, 'Set1');
figure('Position', [440,357,808,441])
l = violin(gammalst, ...
    'mc','','medc','k' , 'facecolor', 'k'); %statecols(4,:));
% l = violin([coefs_all{1}' bias_all{1}'], ...
%     'mc','','medc','k' , 'facecolor', statecols(3,:), 'x', [2,4]);
ylim([-0.05, 0.1])
xlim([0 5])
mymakeaxis('x_label', '', 'y_label', ...
    'Inferred learning rate', 'yticks', -0.05:0.05:0.1, 'xticks', [1,2,3], ...
    'xticklabels', {'Agent 1', 'Agent 2', 'Agent M'}, 'font_size', 20);
legend('off')

%% Aggregate and plot
statecols = brewermap(4, 'Set1');

mean1 = -mean(Ball_arr1, 1);
std1 = std(Ball_arr1, [], 1);
mean2 = -mean(Ball_arr2, 1);
std2 = std(Ball_arr2, [], 1);
meanall = -mean(Ball_arr, 1);
stdall = std(Ball_arr, [], 1);

ymin = -0.2;
ymax = 0.6;
figure;
subplot(131)
errorbar(1:Nlogistic, mean1(Nlogistic:-1:1), std1(1:Nlogistic), 'Color', statecols(1,:), 'LineWidth', 1);
hold on
errorbar(1:Nlogistic, mean2(Nlogistic:-1:1), std2(1:Nlogistic), 'Color', statecols(2,:), 'LineWidth', 1);
errorbar(1:Nlogistic, meanall(Nlogistic:-1:1), stdall(1:Nlogistic), 'Color', statecols(4,:), 'LineWidth', 1);

ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. choice',...
    'xticks', 0:5:10, 'xticklabels', {'0', '-5', '-10'}, 'font_size', 18)



subplot(132)
errorbar(1:Nlogistic, mean1(Nlogistic*2:-1:Nlogistic+1), std1(1:Nlogistic), 'Color', statecols(1,:), 'LineWidth', 1);
hold on
errorbar(1:Nlogistic, mean2(Nlogistic*2:-1:Nlogistic+1), std2(1:Nlogistic), 'Color', statecols(2,:), 'LineWidth', 1);
errorbar(1:Nlogistic, meanall(Nlogistic*2:-1:Nlogistic+1), stdall(1:Nlogistic), 'Color', statecols(4,:), 'LineWidth', 1);

ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. reward',...
    'xticks', 0:5:10, 'xticklabels', {'0', '-5', '-10'}, 'font_size', 18)



subplot(133)
l1 = errorbar(1:Nlogistic, mean1(end:-1:Nlogistic*2+1), std1(1:Nlogistic), 'Color', statecols(1,:), 'LineWidth', 1);
hold on
l2 = errorbar(1:Nlogistic, mean2(end:-1:Nlogistic*2+1), std2(1:Nlogistic), 'Color', statecols(2,:), 'LineWidth', 1);
l3 = errorbar(1:Nlogistic, meanall(end:-1:Nlogistic*2+1), stdall(1:Nlogistic), 'Color', statecols(4,:), 'LineWidth', 1);

ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. choice x reward',...
    'xticks', 0:5:10, 'xticklabels', {'0', '-5', '-10'}, 'font_size', 18)
legend([l1, l2, l3], {'Agent 1', 'Agent 2', 'Agent M'}, 'FontSize', 12);


%% Figure: heatmap of behavior
cols = paperaesthetics;
figure('Position', [440,411,749,387]);
subplot(131)
imagesc(1-block1)
colormap(cols.redbluecolormap2);
ylim([0 1000])
axis xy
mymakeaxis('x_label', 'Trials in block', 'y_label', 'Block #', 'xytitle', 'Agent 1',...
    'xticks', 0:5:15, 'font_size', 18)


subplot(132)
imagesc(1-block2)
colormap(cols.redbluecolormap2);
ylim([0 1000])
axis xy

mymakeaxis('x_label', 'Trials in block', 'y_label', 'Block #', 'xytitle', 'Agent 2',...
    'xticks', 0:5:15, 'font_size', 18)


subplot(133)
imagesc(1-block)
colormap(cols.redbluecolormap2);
ylim([0 1000])
axis xy

mymakeaxis('x_label', 'Trials in block', 'y_label', 'Block #', 'xytitle', 'Agent M',...
    'xticks', 0:5:15, 'font_size', 18)




function block = generate_blocks_with_sigmoid(Nblocks, blocklen, offset, slope, lapse)
% '''
% Generate a session with Nblocks blocks, each block of length
% blocklen trials, governed by transition function with parameters
% (offset, slope, lapse)
% '''
x = 0:blocklen-1;
transfunc = mathfuncs.sigmoid(x, offset, slope, lapse);
block = rand(Nblocks, blocklen) > transfunc;

end



function [choices, outcomes] = unfold_block(block)
outcomes = (1 - block) * 2 - 1;
block = block * 2 - 1;
block(1:2:end,:) = block(1:2:end,:) * -1;
choices = block;

end
