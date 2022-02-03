%% Simulation 1: Simulate sequences of choice and rewards
N = 10000;
penv = 0.8;
pchoice = 0.8;
choices = (rand(N, 1) < pchoice) * 2 - 1;
targets = (rand(N, 1) < penv) * 2 - 1;
outcomes = (choices == targets) * 2 - 1;

% Regression
currchoice = choices(2:end);
prevchoice = choices(1:end-1);
prevrew = outcomes(1:end-1);
pchoicerew = prevchoice .* prevrew;

X = [prevchoice prevrew pchoicerew];
% y = (currchoice > 0) + 1;
% [B, dev, stats] = mnrfit(X, y);
y = currchoice > 0;
b = glmfit(X, y, 'binomial');


%% Simulation 2: a step function from 0.2 to 0.8
rng('shuffle')
N1 = 100;
N2 = 10000;
penv = 0.8;
pchoice1 = 0.2;
pchoice2 = 0.8;

choices1 = (rand(N1, 1) < pchoice1) * 2 - 1;
choices2 = (rand(N2, 1) < pchoice2) * 2 - 1;
choices = [choices1; choices2];
targets = (rand(N1 + N2, 1) < penv) * 2 - 1;
outcomes = (choices == targets) * 2 - 1;

% Regression
currchoice = choices(2:end);
prevchoice = choices(1:end-1);
prevrew = outcomes(1:end-1);
pchoicerew = prevchoice .* prevrew;

X = [prevchoice prevrew pchoicerew];
y = currchoice > 0;
[b,dev,stats] = glmfit(X, y, 'binomial');


%% Simulation 3: sharp switch but alternating blocks
rng('shuffle')
N1 = 100;
N2 = 100;
Nblocks = 20;
penv = 0.8;
pchoice1 = 0.1;
pchoice2 = 0.9;



choicearr = [];
targetarr = [];
outcomearr = [];

for i = 1:Nblocks
    choices1 = (rand(N1, 1) < pchoice1) * 2 - 1;
    choices2 = (rand(N2, 1) < pchoice2) * 2 - 1;
    choices = [choices1; choices2];
    targets = (rand(N1+N2, 1) < penv) * 2 - 1;
    outcomes = (choices == targets) * 2 - 1;
    
    if mod(i, 2) == 0
        choices = choices * -1;
        targets = targets * -1;
    end
    
    choicearr(i,:) = choices;
    targetarr(i,:) = targets;
    outcomearr(i,:) = outcomes;
    
end

choices = reshape(choicearr', [], 1);
targets = reshape(targetarr', [], 1);
outcomes = reshape(outcomearr', [], 1);

% Regression
currchoice = choices(2:end);
prevchoice = choices(1:end-1);
prevrew = outcomes(1:end-1);
pchoicerew = prevchoice .* prevrew;

X = [prevchoice prevrew pchoicerew];
y = currchoice > 0;
[b,dev,stats] = glmfit(X, y, 'binomial');




%% Simulate a sigmoidal transition
N = 10000;
Nswitch = 5000;
Nblocks = 20;
plapse = 0.2;
penv = 0.8;
slope = 0.001;
xvals = (1:N)';
prob = mathfuncs.sigmoid(xvals, Nswitch, slope, plapse);
choices = (rand(N, Nblocks) < prob) * 2 - 1;

targets = (rand(N, 1) < penv) * 2 - 1;
outcomes = (choices == targets) * 2 - 1;

% Regression
currchoice = choices(2:end);
prevchoice = choices(1:end-1);
prevrew = outcomes(1:end-1);
pchoicerew = prevchoice .* prevrew;

X = [prevchoice prevrew pchoicerew];
y = currchoice > 0;
[b,dev,stats] = glmfit(X, y, 'binomial');



%% Simulation 4: simulate multiple sigmoidal transitions
opts.N = 20;
opts.Nswitch = 5;
opts.Nblocks = 2000;
opts.plapse = 0.2;
opts.penv = 0.8;
opts.slope = 2;
opts.Nreps = 100;
nswitchlst = 2:18;
slopelst = linspace(0.01, 5, 10);
b_all = [];
b2arr = [];
b4arr = [];
for j = 1:numel(nswitchlst)
    for k = 1:numel(slopelst)
        opts.Nswitch = nswitchlst(j);
        opts.slope = slopelst(k);
        for i = 1:opts.Nreps
            b = perform_simulation(opts);
            b_all(i,:) = b;
        end
        fprintf('Nswitch = %d, slope = %.2f: Mean b2 = %.2f, mean b4 = %.2f\n', opts.Nswitch, opts.slope, mean(b_all(:,2)), mean(b_all(:,4)));
        b2arr(j,k) = mean(b_all(:,2));
        b4arr(j,k) = mean(b_all(:,4));
    end
end

%%
figure;
subplot(121)
imagesc(b2arr)
colorbar

subplot(122)
imagesc(b4arr)
colorbar


function b = perform_simulation(opts)
N = opts.N;
Nswitch = opts.Nswitch;
Nblocks = opts.Nblocks;
plapse = opts.plapse;
penv = opts.penv;
slope = opts.slope;
xvals = (1:N)';
prob = mathfuncs.sigmoid(xvals, Nswitch, slope, plapse);
choicearr = (rand(N, Nblocks) < prob) * 2 - 1;
targetarr = (rand(N, 1) < penv) * 2 - 1;
outcomearr = (choicearr == targetarr) * 2 - 1;

%flip the even trials
choicearr(:,1:2:end) = choicearr(:,1:2:end) * -1;
targetarr(:,1:2:end) = targetarr(:,1:2:end) * -1;

choices = reshape(choicearr, [], 1);
targets = reshape(targetarr, [], 1);
outcomes = reshape(outcomearr, [], 1);

% Regression
currchoice = choices(2:end);
prevchoice = choices(1:end-1);
prevrew = outcomes(1:end-1);
pchoicerew = prevchoice .* prevrew;

X = [prevchoice prevrew pchoicerew];
y = currchoice > 0;
[b,dev,stats] = glmfit(X, y, 'binomial');

end


