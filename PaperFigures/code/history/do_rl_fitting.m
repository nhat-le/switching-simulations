function res = do_rl_fitting(choices, outcomes)
    choices_flat = reshape(choices', [], 1);
    outcomes_flat = reshape(outcomes', [], 1);
    outcomes_flat = (outcomes_flat + 1) / 2;
    mdl = rl.fit(choices_flat, outcomes_flat);
    res = mdl.params;
end