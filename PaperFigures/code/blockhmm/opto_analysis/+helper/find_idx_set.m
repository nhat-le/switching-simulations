function [idxperiod, idxpower] = find_idx_set(modecounts, power_criterion, ...
    period_criterion)

if strcmp(power_criterion, 'none')
    idxpower = [1,2];
else
    idxpower = find(contains(modecounts(1).power_criteria, power_criterion));
end

if strcmp(period_criterion, 'none')
    idxperiod = [1,2];
else
    idxperiod = find(contains(modecounts(1).period_criteria, period_criterion));
end

end