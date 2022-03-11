%% Data generation
seed = 124;
rng(seed);
Niter = 10;
Nlogistic = 10;
Nsubblock = 50;

h = waitbar(0);

gammalst = [];
biaslst = [];
Ball_arr = [];
Ball_arr1 = [];
Ball_arr2 = [];

for i = 1:Niter
    % Load the data
    waitbar(i/Niter, h)
    % Qlearning data
    load(sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/pitfalls/qlearning_gamma_0.1_eps0.1_iter%d.mat', ...
        i-1));
    choiceslst1 = choices;
    outcomeslst1 = outcomes;
    choicearr1 = choicearr;
    outcomearr1 = outcomearr;
    
      
    % Inf-based data
    load(sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata/pitfalls/inf-based_psw_0.45_pr0.99_iter%d.mat',...
        i-1));
    choiceslst2 = choices;
    outcomeslst2 = outcomes;
    choicearr2 = choicearr;
    outcomearr2 = outcomearr;

    
    % Mixed data
    choicearrM = mix_arrays(choicearr1, choicearr2, Nsubblock);
    outcomearrM = mix_arrays(outcomearr1, outcomearr2, Nsubblock);
    

    %RL fitting
    res = do_rl_fitting(choicearrM, outcomearrM);
    res1 = do_rl_fitting(choicearr1, outcomearr1);
    res2 = do_rl_fitting(choicearr2, outcomearr2);
    
    gammalst(end+1,:) = [res1(1), res2(1), res(1)];
    biaslst(end+1,:) = [res1(3), res2(3), res(3)];
    
    % Logistic regression
%     B_all = do_logistic_regression(block, Nlogistic);
%     B_all1 = do_logistic_regression(block1, Nlogistic);
%     B_all2 = do_logistic_regression(block2, Nlogistic);
%     
%     Ball_arr(end+1,:) = B_all;
%     Ball_arr1(end+1,:) = B_all1;
%     Ball_arr2(end+1,:) = B_all2;

end
close(h)


%% Figure: RL summary
statecols = brewermap(4, 'Set1');
figure('Position', [440,357,808,441])
l = violin(gammalst, ...
    'mc','','medc','k' , 'facecolor', 'k'); %statecols(4,:));
% l = violin([coefs_all{1}' bias_all{1}'], ...
%     'mc','','medc','k' , 'facecolor', statecols(3,:), 'x', [2,4]);
% ylim([-0.05, 0.1])
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
errorbar(1:Nlogistic, mean1(Nlogistic:-1:1), std1(Nlogistic:-1:1), 'Color', statecols(1,:), 'LineWidth', 1);
hold on
errorbar(1:Nlogistic, mean2(Nlogistic:-1:1), std2(Nlogistic:-1:1), 'Color', statecols(2,:), 'LineWidth', 1);
errorbar(1:Nlogistic, meanall(Nlogistic:-1:1), stdall(Nlogistic:-1:1), 'Color', statecols(4,:), 'LineWidth', 1);

ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. choice',...
    'xticks', 0:5:10, 'xticklabels', {'0', '-5', '-10'}, 'font_size', 18)



subplot(132)
errorbar(1:Nlogistic, mean1(Nlogistic*2:-1:Nlogistic+1), std1(Nlogistic*2:-1:Nlogistic+1), 'Color', statecols(1,:), 'LineWidth', 1);
hold on
errorbar(1:Nlogistic, mean2(Nlogistic*2:-1:Nlogistic+1), std2(Nlogistic*2:-1:Nlogistic+1), 'Color', statecols(2,:), 'LineWidth', 1);
errorbar(1:Nlogistic, meanall(Nlogistic*2:-1:Nlogistic+1), stdall(Nlogistic*2:-1:Nlogistic+1), 'Color', statecols(4,:), 'LineWidth', 1);

ylim([ymin ymax])
mymakeaxis('x_label', 'Trials', 'y_label', 'Coefficient', 'xytitle', 'Prev. reward',...
    'xticks', 0:5:10, 'xticklabels', {'0', '-5', '-10'}, 'font_size', 18)



subplot(133)
l1 = errorbar(1:Nlogistic, mean1(end:-1:Nlogistic*2+1), std1(end:-1:Nlogistic*2+1), 'Color', statecols(1,:), 'LineWidth', 1);
hold on
l2 = errorbar(1:Nlogistic, mean2(end:-1:Nlogistic*2+1), std2(end:-1:Nlogistic*2+1), 'Color', statecols(2,:), 'LineWidth', 1);
l3 = errorbar(1:Nlogistic, meanall(end:-1:Nlogistic*2+1), stdall(end:-1:Nlogistic*2+1), 'Color', statecols(4,:), 'LineWidth', 1);

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



function ll = neg_choice_likelihood(response, feedback, p)
% beta: sharpness parameter of the sigmoid
alpha = p(1);
beta = p(2);
offset = p(3);

[v0, v1] = rl.create_value_arr(alpha, response, feedback);
prob = sigmoid(v1 - v0, beta, offset);
if sum(size(response) == size(prob)) == 0
    prob = prob';
end
responsebinary = response > 0;
ll = responsebinary .* log(prob) + (1-responsebinary) .* log(1 - prob);
ll = -sum(ll);


end


function y = sigmoid(x, beta, offset)
y = 1 ./ (1 + exp(-(x - offset) * beta));

end


function mixedarr = mix_arrays(arr1, arr2, Nsubblock)
% arr1 and arr2 should have the same size
% Nblocks x T
% Create another array with Nblocks x T by mixing arr1 and arr2
% in sub-blocks of lengths given by Nsubblock

[Nblocks1, T1] = size(arr1);
[Nblocks2, T2] = size(arr2);

assert(Nblocks1 == Nblocks2);
assert(T1 == T2);

mixedarr = [];

for i = 1:floor(Nblocks1 / Nsubblock / 2)
    idxstart = (i - 1) * Nsubblock + 1;
    idxend = i * Nsubblock;
    mixedarr = [mixedarr; arr1(idxstart:idxend, :)];
    mixedarr = [mixedarr; arr2(idxstart:idxend, :)];

end


end
