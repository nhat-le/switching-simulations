function res = repeated_rl_fitting(choices, outcomes)
% choices, outcomes: size Nblocks x Ntrialsperblock
% Performs resampling and fitting of the RL model
%
choices(1:2:end) = choices(1:2:end) * -1;
% Resample rows of choices and outcomes
assert(size(choices,1) == size(outcomes,1));
assert(size(choices,2) == size(outcomes,2));

idx = randsample(1:size(choices,1), size(choices,1), true);
choice_resampled = choices(idx,:);
choice_resampled(1:2:end) = choice_resampled(1:2:end) * -1;
outcomes_resampled = outcomes(idx,:);

choices_flat = reshape(choice_resampled', [], 1);
outcomes_flat = reshape(outcomes_resampled', [], 1);
outcomes_flat = (outcomes_flat + 1) / 2;
mdl = rl.fit(choices_flat, outcomes_flat);
res = mdl.params;

end