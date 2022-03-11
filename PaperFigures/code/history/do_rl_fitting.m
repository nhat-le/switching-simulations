function res = do_rl_fitting(choices, outcomes)
    % choices: array of size Nblocks x T, values 0 or 1
    % outcomes: array of size Nblocks x T, values 0 or 1
    assert(max(choices(:)) == 1)
    assert(max(outcomes(:)) == 1)
    assert(min(choices(:)) == 0)
    assert(min(outcomes(:)) == 0)
    
    
    % rl.fit works with +1/-1 choices, so do the transformation here
    choices = choices * 2 - 1;
    
    [Nblocks1, T1] = size(choices);
    [Nblocks2, T2] = size(outcomes);
    assert(Nblocks1 == Nblocks2);
    assert(T1 == T2);
  
    choices_flat = reshape(choices', [], 1);
    outcomes_flat = reshape(outcomes', [], 1);
    
    filt_idx = (~isnan(choices_flat) & ~isnan(outcomes_flat));
    choices_flat = choices_flat(filt_idx);
    outcomes_flat = outcomes_flat(filt_idx);
    
    assert(sum(isnan(choices_flat)) == 0);
    assert(sum(isnan(outcomes_flat)) == 0);
    
    
    mdl = rl.fit(choices_flat, outcomes_flat);
    res = mdl.params;
end