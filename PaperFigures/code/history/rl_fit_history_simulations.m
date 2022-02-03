% load('shuffle_simulation_results_v3.mat');
% 
% rng('shuffle')
% Ntrials = size(block, 1) * size(block, 2);
% [Nrows, Ncols] = size(Xmat);
% choices1 = double(Xmat(:,1));
% outcomes1 = double((Xmat(:,Ncols/3+1) + 1) / 2);
% mdl1 = rl.fit(choices1, outcomes1);
% 
% choices2 = double(Xmat2(:,1));
% outcomes2 = double((Xmat2(:,Ncols/3+1) + 1) / 2);
% mdl2 = rl.fit(choices2, outcomes2);
% 
% 
% choices3 = double(Xmat3(:,1));
% outcomes3 = double((Xmat3(:,Ncols/3+1) + 1) / 2);
% mdl3 = rl.fit(choices3, outcomes3);
% 
% fprintf('--- Results: ---\n');
% disp(mdl1.params)
% disp(mdl2.params)
% disp(mdl3.params)

%% Data generation
Nblocks = 1000;
blocklen = 20;

offset1 = 5;
offset2 = 5;
slope1 = 0.6;
slope2 = 2;
lapse1 = 0;
lapse2 = 0;

offset3 = 4.8273;
slope3 = 0.8214;
lapse3 = 0;
block1 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset1, slope1, lapse1);
block2 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset2, slope2, lapse2);
block = [block1; block2];


shuffleblock = generate_blocks_with_sigmoid(2 * Nblocks, blocklen, offset3, slope3, lapse3);
 
%% Repeated RL fitting
Ntrials = 100;
coefslst1 = repeated_rl_fitting_coefs(block, Ntrials);
coefslst2 = repeated_rl_fitting_coefs(shuffleblock, Ntrials);
% coefslst3 = repeated_rl_fitting_coefs(shuffleblock, Ntrials);

%%
mean1 = mean(coefslst1);
mean2 = mean(coefslst2);
std1 = std(coefslst1', [], 1);
std2 = std(coefslst2', [], 1);
figure;
subplot(121)
plot(coefslst1)
hold on
plot(coefslst2)

subplot(122)
errorbar(1:2, [mean1, mean2], [std1, std2], 'o-')
xlim([0, 3])


%% Logistic regression with error bars
Ntrials = 10;
N = 15;
B_all = do_logistic_regression_with_error_bars(block, Ntrials, N);

%%
B_all2 = do_logistic_regression_with_error_bars(shuffleblock, Ntrials, N);

%% Aggregate and plot
mean1 = mean(B_all, 1);
std1 = std(B_all, [], 1);
mean2 = mean(B_all2, 1);
std2 = std(B_all2, [], 1);


figure;
subplot(131)
errorbar(1:N, mean1(1:N), std1(1:N))
hold on
errorbar(1:N, mean2(1:N), std2(1:N))

ylim([-1 1])

subplot(132)
errorbar(1:N, mean1(N+1:2*N), std1(1:N))
hold on
errorbar(1:N, mean2(N+1:2*N), std2(1:N))
ylim([-1 1])


subplot(133)
errorbar(1:N, mean1(2*N+1:end), std1(1:N))
hold on
errorbar(1:N, mean2(2*N+1:end), std2(1:N))
ylim([-1 1])

function B_all = do_logistic_regression_with_error_bars(block, Ntrials, N)
[choicelst1, outcomelst1] = unfold_block(block);
B_all = [];
h = waitbar(0);
for i = 1:Ntrials
    waitbar(i/Ntrials, h);
    [Xmat1, y1] = make_Xy(choicelst1, outcomelst1, N, true);

    % logistic regression
    B = mnrfit(Xmat1, (y1 > 0) + 1);
    B = B(2:end); %first coef is intercept
    B_all(end+1,:) = B;
end
close(h);
end

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



function coefslst = repeated_rl_fitting_coefs(blockarr, Ntrials)
[choiceslst, outcomeslst] = unfold_block(blockarr);
coefslst = [];
h = waitbar(0);
for i = 1:Ntrials
    waitbar(i/Ntrials, h);
    res = repeated_rl_fitting(choiceslst, outcomeslst);
    coefslst(end+1) = res(1);
end
close(h);

end



function [choices, outcomes] = unfold_block(block)
outcomes = (1 - block) * 2 - 1;
block = block * 2 - 1;
block(1:2:end,:) = block(1:2:end,:) * -1;
choices = block;

end
