function [choices, outcomes] = unfold_block(block)
% block: size Nblocks x T, 0 or 1
% returns: choices: size Nblocks x T, -1 or 1
% outcomes: size Nblocks x T, -1 or 1

outcomes = (1 - block) * 2 - 1;
block = block * 2 - 1;
block(1:2:end,:) = block(1:2:end,:) * -1;
choices = block;

end