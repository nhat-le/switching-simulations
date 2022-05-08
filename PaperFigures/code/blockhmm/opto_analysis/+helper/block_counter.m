function counts = block_counter(arr, idxperiod, idxpower)
% arr: a cell with n power x n area x n periods
% idxperiod: array, subset of periods of interest
% idxpower: array, subset of powers of interest
% counts: nareas x nstates array, with the counts of n blocks in each state
% per area
[npower, nareas, nperiods] = size(arr);
assert(max(idxperiod) <= nperiods);
assert(max(idxpower) <= npower);

nstates = numel(arr{1,1,1});

counts = zeros(nareas, nstates);

for i = 1:numel(idxpower)
    poweridx = idxpower(i);
    for j = 1:nareas
        for k = 1:numel(idxperiod)
            periodidx = idxperiod(k);
            count_condition = arr{poweridx, j, periodidx};
            
            counts(j, :) = counts(j, :) + count_condition;
            
        end
    end
end

end
