function [Xmat, y] = make_Xy_flat(choices, outcomes, N)
% Inputs: % choicearr: array of choices, size Nblocks x T, -1 or 1
% outcomes: array of outcomes, size Nblocks x T, -1 or 1
% N: how many trials back to go for logistic regression
% Prepares X, y data ready for logistic regression
% Xmat has dimension (Ntrials - N) x 3N

assert(max(choices(:)) == 1);
assert(min(choices(:)) == -1);
assert(max(outcomes(:)) == 1);
assert(min(outcomes(:)) == -1);


if isempty(N)
    N = 1;
end

choices_flat = reshape(choices', [], 1);
outcomes_flat = reshape(outcomes', [], 1);

% Eliminate nan's
filt_idx = (~isnan(choices_flat) & ~isnan(outcomes_flat));
choices_flat = choices_flat(filt_idx);
outcomes_flat = outcomes_flat(filt_idx);

assert(sum(isnan(choices_flat)) == 0)
assert(sum(isnan(outcomes_flat)) == 0);

outcomesXchoices = choices_flat .* outcomes_flat;

assert(max(choices_flat(:)) == 1);
assert(min(choices_flat(:)) == -1);
assert(max(outcomes_flat(:)) == 1);
assert(min(outcomes_flat(:)) == -1);
assert(max(outcomesXchoices(:)) == 1);
assert(min(outcomesXchoices(:)) == -1);


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

assert(size(Xmat, 2) == 3 * N);
assert(rank(Xmat) == 3 * N); % ensures convergence of the algorithm


end