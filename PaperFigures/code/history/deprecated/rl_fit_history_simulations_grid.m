%% Data generation
Nblocks = 1000;
blocklen = 15;

offset1 = 5;
offset2 = 4;
slope1 = 0.6;
slope2 = 2;
lapse1 = 0.1;
lapse2 = 0.1;

seed = 124;
rng(seed);
offsetlst = ones(1,1) * 4.6;
slopelst = ones(1,1) * 0.85;
% offsetlst = ones(1,10) * 4.6;
% slopelst = ones(1,10) * 0.85;
% offsetlst = rand(1, 25)*1 + 4.5; %linspace(4.5, 5.5, 5);
% slopelst = rand(1, 25) * 0.3 + 0.75; %linspace(0.7, 0.9, 5);

lapse3 = 0.1;
block1 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset1, slope1, lapse1);
block2 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset2, slope2, lapse2);

interleave = 1;
block = [block1; block2];
Nsubblock = 250;
idxall = [];
if interleave
    % Interleave blocks 1 and 2 in subblocks of 200
    idx = [1:Nsubblock,  (Nblocks + 1) : (Nblocks + Nsubblock)];
    for i = 1:(Nblocks / Nsubblock)
        idxall = [idxall idx];
        idx = idx + Nsubblock;
        
    end 
    block = block(idxall, :);
end

coefs_all = {};
bias_all = {};
Niter = 100;

for i = 1:numel(offsetlst)
    fprintf('Fitting element %d\n', i);
%     shuffleblock = generate_blocks_with_sigmoid(2 * Nblocks, blocklen, offsetlst(i), slopelst(i), lapse3);
%     shuffleblocks{i} = shuffleblock;
    [coefslst, biaslst] = repeated_rl_sim_and_fitting(2*Nblocks, blocklen, offsetlst(i), slopelst(i), lapse3, Niter);
%     coefslst = repeated_rl_fitting_coefs(shuffleblock, Ntrials);
    coefs_all{i} = coefslst;
    bias_all{i} = biaslst;
end

% Repeated RL fitting
fprintf('Fitting the mixed behavior..\n')
[coefslst1, biaslst1] = repeated_rl_sim_and_fitting_mixed(Nblocks, blocklen, offset1, slope1, lapse1, ...
    offset2, slope2, lapse2, Niter, Nsubblock, interleave);


%%
arr = cell2mat(coefs_all');
meanall = mean(arr, 2);
stdall = std(arr, [], 2);

diffmean = (meanall - mean(coefslst1)).^2;
id = argmin(diffmean);
fprintf('Best id = %d, offset = %.2f, slope = %.2f\n', id, offsetlst(id), slopelst(id));
diffmean(id) = inf;
id2 = argmin(diffmean);
fprintf('Second best id = %d, offset = %.2f, slope = %.2f\n', id2, offsetlst(id2), slopelst(id2));


figure;
subplot(121)
errorbar(1:size(arr, 1), meanall, stdall, 'o-');
hold on
errorbar(int16(size(arr, 1) / 2), mean(coefslst1), std(coefslst1, [], 2), 'o-')

subplot(122)
plot(1, coefslst1, 'ro')
hold on
plot(2, coefs_all{1}, 'bo')

%% save
blockalt = generate_blocks_with_sigmoid(2*Nblocks, blocklen, offsetlst(1), slopelst(1), lapse3);

% save('rlmodel_sim_020122_v5.mat', 'seed', 'Nblocks', 'blocklen', 'Niter',...
%     'offset1', 'offset2', 'slope1', 'slope2', 'lapse1', 'lapse2', 'lapse3', ...
%     'slopelst', 'offsetlst', 'coefs_all', 'coefslst1', 'block', 'blockalt',...
%     'biaslst1', 'bias_all');



%% Logistic regression with error bars
% Ntrials = 10;
% N = 15;
% B_all = do_logistic_regression_with_error_bars(block, Ntrials, N);
% 
% %%
% B_all2 = do_logistic_regression_with_error_bars(shuffleblock, Ntrials, N);
% 

%%
seed = 123;
rng(seed);
N = 15;
Nblocks = 1000;
blocklen = 15;
Niter = 10;

offset1 = 5;
offset2 = 4;
slope1 = 0.6;
slope2 = 2;
lapse1 = 0.1;
lapse2 = 0.1;
lapse3 = 0.1;

Nsubblock = 200;
interleave = 1;

offsetlst = ones(1,1) * 4.6;
slopelst = ones(1,1) * 0.84;
% offsetlst = rand(1, 25)*2 + 4; %linspace(4.5, 5.5, 5);
% slopelst = rand(1, 25) * 0.5 + 0.5; %linspace(0.7, 0.9, 5);


fprintf('Fitting the mixed behavior ...\n');
B_all1 = repeated_logistic_sim_and_fitting_mixed(Nblocks, N, blocklen, offset1, slope1, lapse1, ...
    offset2, slope2, lapse2, Niter, Nsubblock, interleave);

B_all_cell = {};
for i = 1:numel(offsetlst)
    fprintf('Fitting step %d of %d...\n', i, numel(offsetlst))
    B_all2 = repeated_logistic_sim_and_fitting(2*Nblocks, N, blocklen, offsetlst(i), slopelst(i), lapse3, Niter);
    B_all_cell{i} = B_all2;
end

%%

Bmean1 = mean(B_all1, 1);
Bdiffs = [];
for i = 1:numel(B_all_cell)
    B_all2 = B_all_cell{i};
    Bmean = mean(B_all2, 1);
    
    Bdiff = sum((Bmean - Bmean1).^2);
    Bdiffs(i) = Bdiff;
end

id = argmin(Bdiffs);
fprintf('Best id = %d, offset = %.2f, slope = %.2f\n', id, offsetlst(id), slopelst(id));
Bdiffs(id) = inf;
id2 = argmin(Bdiffs);
fprintf('Second best id = %d, offset = %.2f, slope = %.2f\n', id2, offsetlst(id2), slopelst(id2));



%% Aggregate and plot
for i = 1
    mean1 = mean(B_all1, 1);
    std1 = std(B_all1, [], 1);
    B_all2 = B_all_cell{i};
    mean2 = mean(B_all2, 1);
    std2 = std(B_all2, [], 1);


    figure;
    subplot(131)
    errorbar(1:N, mean1(1:N), std1(1:N))
    hold on
    errorbar(1:N, mean2(1:N), std2(1:N))

    ylim([-3 1])

    subplot(132)
    errorbar(1:N, mean1(N+1:2*N), std1(1:N))
    hold on
    errorbar(1:N, mean2(N+1:2*N), std2(1:N))
    ylim([-3 1])


    subplot(133)
    errorbar(1:N, mean1(2*N+1:end), std1(1:N))
    hold on
    errorbar(1:N, mean2(2*N+1:end), std2(1:N))
    ylim([-3 1])
end

%% Save processed data
% save('logistic_sim_013022_v1.mat',  'seed', 'N', 'Nblocks', 'blocklen', 'Niter',...
%     'offset1', 'offset2', 'slope1', 'slope2', 'lapse1', 'lapse2', 'lapse3', ...
%     'slopelst', 'offsetlst', 'B_all1', 'B_all_cell');





function B_all = repeated_logistic_sim_and_fitting(Nblocks, N, ...
    blocklen, offset, slope, lapse, Niter)

h = waitbar(0);
B_all = [];
for i = 1:Niter
    waitbar(i/Niter, h);
    block = generate_blocks_with_sigmoid(Nblocks, blocklen, offset, slope, lapse);
    B = do_logistic_regression(block, N);
    B_all(end+1,:) = B;    
end
close(h);

end


function B_all = repeated_logistic_sim_and_fitting_mixed(Nblocks, N, ...
    blocklen, offset1, slope1, lapse1, offset2, slope2, lapse2, Niter, Nsubblock, interleave)

h = waitbar(0);
B_all = [];
for i = 1:Niter
    waitbar(i/Niter, h);
    block1 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset1, slope1, lapse1);
    block2 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset2, slope2, lapse2);
    block = [block1; block2];
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
    B = do_logistic_regression(block, N);
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


