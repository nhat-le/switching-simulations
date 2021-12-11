% A script for comparing the params between different conditions
IBfiles = dir('EGreedyInferenceBasedAgent*.mat');
Qfiles = dir('EGreedyQLearning*.mat');

%% Collect the 'best' param in terms of foraging efficiency for IB
prmaxIB = [];
psmaxIB = [];
rlowsIB = [];
for i = 1:5
    load(fullfile(IBfiles(i).folder, IBfiles(i).name));
    rlowsIB(i) = rlow;
    [idx, idy] = argmaxArray(efflist);
    prmaxIB(i) = prewlst(idy);
    psmaxIB(i) = pswitchlst(idx);
    
    figure;
    imagesc(efflist);
    hold on
    plot(idy, idx, 'x')
    axis xy
    
end

%% Collect the 'best' param in terms of foraging efficiency for Q
epsmaxQ = [];
gammamaxQ = [];
rlowsQ = [];
for i = 1:5
    load(fullfile(Qfiles(i).folder, Qfiles(i).name));
    rlowsQ(i) = rlow;
    [idx, idy] = argmaxArray(efflist);
    epsmaxQ(i) = epslst(idy);
    gammamaxQ(i) = gammalst(idx);
    figure;
    imagesc(efflist);
    hold on
    plot(idy, idx, 'x')
    axis xy
    
end


%%
figure;
scatter(prmaxIB, psmaxIB, 30, rlowsIB + 0.1, 'o');


%%
figure;
scatter(epsmaxQ, gammamaxQ, 30, rlowsQ + 0.1, 'o');
