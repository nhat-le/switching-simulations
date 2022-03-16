function [coefslst, biaslst] = repeated_rl_sim_and_fitting_mixed(Nblocks, blocklen, offset1, slope1, lapse1, ...
    offset2, slope2, lapse2, Niter, Nsubblock, interleave)
% choices, outcomes: size Nblocks x Ntrialsperblock
% Performs resampling and fitting of the RL model
%
coefslst = [];
biaslst = [];
h = waitbar(0);
for i = 1:Niter
    waitbar(i/Niter, h);
    block1 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset1, slope1, lapse1);
    block2 = generate_blocks_with_sigmoid(Nblocks, blocklen, offset2, slope2, lapse2);
    block = [block1; block2];
    
    idxall = [];
    if interleave
        % Interleave blocks 1 and 2 in subblocks of 200
        idx = [1:Nsubblock,  (Nblocks + 1) : (Nblocks + Nsubblock)];
        for j = 1:(Nblocks / Nsubblock)
            idxall = [idxall idx];
            idx = idx + Nsubblock;

        end 
        block = block(idxall, :);
    end

    [choiceslst, outcomeslst] = unfold_block(block);
    res = do_rl_fitting(choiceslst, outcomeslst);
    coefslst(end+1) = res(1);
    biaslst(end+1) = res(3);
end
close(h);

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