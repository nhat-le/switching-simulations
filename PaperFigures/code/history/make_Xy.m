function [Xmat, y] = make_Xy(choices, outcomes, N, resample)
% Inputs: choices: array of size Nblocks x T, -1 or 1
% outcomes: array of size Nblocks x T, -1 or 1
% N: how many trials back to go for logistic regression
% resample: if True, will resample from the arrays

assert(max(choices(:)) == 1);
assert(min(choices(:)) == -1);
assert(max(outcomes(:)) == 1);
assert(min(outcomes(:)) == -1);


if isempty(N)
    N = 1;
end

% Resample if requested
if resample
    choices(1:2:end) = choices(1:2:end) * -1;
    % Resample rows of choices and outcomes
    assert(size(choices,1) == size(outcomes,1));
    assert(size(choices,2) == size(outcomes,2));

    idx = randsample(1:size(choices,1), size(choices,1), true);
    choice_resampled = choices(idx,:);
    choice_resampled(1:2:end) = choice_resampled(1:2:end) * -1;
    outcomes_resampled = outcomes(idx,:); 
else
    choice_resampled = choices;
    outcomes_resampled = outcomes;
end

choices_flat = reshape(choice_resampled', [], 1);
outcomes_flat = reshape(outcomes_resampled', [], 1);
choices_flat = choices_flat(~isnan(choices_flat));
outcomes_flat = outcomes_flat(~isnan(outcomes_flat));
outcomesXchoices = choices_flat .* outcomes_flat;

y = choices_flat(N + 1 : end);
choicehistory = [];
outcomehistory = [];
choiceoutcomehistory = [];

for i = 1:N
    X1 = choices_flat(i: end-N+i-1);
    X2 = outcomes_flat(i: end-N+i-1);
    X3 = outcomesXchoices(i: end-N+i-1);
    choicehistory(end + 1,:) = X1;
    outcomehistory(end + 1,:) = X2;
    choiceoutcomehistory(end + 1,:) = X3;
    
end

Xmat = [choicehistory; outcomehistory; choiceoutcomehistory]';


end