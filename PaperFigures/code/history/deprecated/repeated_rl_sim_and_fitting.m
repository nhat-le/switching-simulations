function [coefslst, biaslst] = repeated_rl_sim_and_fitting(Nblocks, blocklen, offset, slope, lapse, Niter)
% choices, outcomes: size Nblocks x Ntrialsperblock
% Performs resampling and fitting of the RL model
%
coefslst = [];
biaslst = [];
h = waitbar(0);
for i = 1:Niter
    waitbar(i/Niter, h);
    block = generate_blocks_with_sigmoid(Nblocks, blocklen, offset, slope, lapse);
    [choiceslst, outcomeslst] = unfold_block(block);
    res = do_rl_fitting(choiceslst, outcomeslst);
    coefslst(end+1) = res(1);
    biaslst(end+1) = res(3);
end
close(h);

end

function res = do_rl_fitting(choices, outcomes)
    choices_flat = reshape(choices', [], 1);
    outcomes_flat = reshape(outcomes', [], 1);
    outcomes_flat = (outcomes_flat + 1) / 2;
    mdl = rl.fit(choices_flat, outcomes_flat);
    res = mdl.params;
end


function [choices, outcomes] = unfold_block(block)
outcomes = (1 - block) * 2 - 1;
block = block * 2 - 1;
block(1:2:end,:) = block(1:2:end,:) * -1;
choices = block;

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