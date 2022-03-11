load('rl_shuffle_simulation_results_v1.mat');

rng('shuffle')
% Ntrials = size(block, 1) * size(block, 2);
[Nrows, Ncols] = size(Xmat1);
% choices1 = double(Xmat1(:,1));
% outcomes1 = double((Xmat1(:,Ncols/3+1) + 1) / 2);
choices1a = reshape(choicelstA', [], 1);
outcomes1a = reshape(outcomelstA', [], 1);
outcomes1a = (outcomes1a + 1) / 2;
mdl1 = rl.fit(choices1a, outcomes1a);


choices2 = double(Xmat2(:,1));
outcomes2 = double((Xmat2(:,Ncols/3+1) + 1) / 2);
mdl2 = rl.fit(choices2, outcomes2);

fprintf('--- Results: ---\n');
disp(mdl1.params)

disp(mdl2.params)

%%
coefslst = [];
h = waitbar(0);
Ntrials = 10;
for i = 1:Ntrials
    waitbar(i/Ntrials, h);
    res = repeated_rl_fitting(choicelstA, outcomelstA);
    coefslst(end+1) = res(1);
end
close(h);

%%
coefslst2 = [];
h = waitbar(0);
for i = 1:Ntrials
    waitbar(i/Ntrials, h);
    res = repeated_rl_fitting(choicelst2, outcomelst2);
    coefslst2(end+1) = res(1);
end
close(h);

%%
figure;
plot(coefslst)
hold on
plot(coefslst2)

%% Logistic regression with error bars
N = 15;

B_all = [];
h = waitbar(0);
for i = 1:Ntrials
    waitbar(i/Ntrials, h);
    [Xmat1, y1] = make_Xy(choicelstA, outcomelstA, N, true);

    % logistic regression
    B = mnrfit(Xmat1, (y1 > 0) + 1);
    B = B(2:end); %first coef is intercept
    B_all(end+1,:) = B;
end
close(h);

%%
B_all2 = [];
h = waitbar(0);
for i = 1:Ntrials
    waitbar(i/Ntrials, h);
    [Xmat2, y2] = make_Xy(choicelst2, outcomelst2, N, true);
    
    % logistic regression
    B = mnrfit(Xmat2, (y1 > 0) + 1);
    B = B(2:end); %first coef is intercept
    B_all2(end+1,:) = B;
end
close(h);

%% Aggregate and plot
mean1 = mean(B_all, 1);
std1 = std(B_all, [], 1);
errorbar(1:size(B_all, 2), mean1, std1)
hold on

mean2 = mean(B_all2, 1);
std2 = std(B_all2, [], 1);
errorbar(1:size(B_all2, 2), mean2, std2)


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





