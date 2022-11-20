function out = get_flags_300um(animalinfo, animalID)
% modified from: helper.get_flags
% include criterion for 300um flag (looking for 'frontal300um' instead)

% animalinfo: struct saved by previous analysis
% animalid: id from 1 - n of animal of interest
% returns: out, structure containing all flags extracted from the animalinfo
% power_flags = {low_power_flag, high_power_flag};
% area_flags = {frontal_flag, motor_flag, visual_flag, rsc_flag};
% period_flags = {choice_flag, outcome_flag};



% Parse information for the animal for each session
opto_flag = cellfun(@(x) sum(x) > 0, animalinfo(animalID).opto);

% no power information if not f26, 27, 29, 32
if ~ismember(lower(animalinfo(animalID).animal), {'f26', 'f27', 'f29', 'f32'})
    low_power_flag = ones(1, numel(animalinfo(animalID).power));
    high_power_flag = 1 - low_power_flag;
else
    low_power_flag = cellfun(@(x) strcmp(x, 'Low'), animalinfo(animalID).power);
    high_power_flag = cellfun(@(x) strcmp(x, 'High'), animalinfo(animalID).power);
end
frontal_flag = cellfun(@(x) contains(x, 'Frontal300') & ~contains(x, 'Motor300') & ~contains(x, 'bad'), ...
    animalinfo(animalID).areas);
motor_flag = cellfun(@(x) contains(x, 'Motor300') & ~contains(x, 'Visual300') & ~contains(x, 'bad'), ...
    animalinfo(animalID).areas);
rsc_flag = cellfun(@(x) contains(x, 'Rsc300') & ~contains(x, 'Motor300') & ~contains(x, 'bad'), ...
    animalinfo(animalID).areas);
visual_flag = cellfun(@(x) contains(x, 'Visual300') & ~contains(x, 'Motor300') & ~contains(x, 'bad'), ...
    animalinfo(animalID).areas);
choice_flag = cellfun(@(x) contains(x, 'Choice'), ...
    animalinfo(animalID).period);
outcome_flag = cellfun(@(x) contains(x, 'Outcome'), ...
    animalinfo(animalID).period);

% make sure all flags are of the same size
assert(max([numel(opto_flag), numel(low_power_flag), numel(frontal_flag),...
    numel(motor_flag), numel(rsc_flag), numel(visual_flag)]) == numel(rsc_flag));
assert(min([numel(opto_flag), numel(low_power_flag), numel(frontal_flag),...
    numel(motor_flag), numel(rsc_flag), numel(visual_flag)]) == numel(rsc_flag));
% assert(sum(low_power_flag + high_power_flag ~= opto_flag) == 0);
% Now we will gather all blocks with the right opto information

% other checks for consistency
assert(max(low_power_flag + high_power_flag) == 1);
assert(max(frontal_flag + visual_flag + rsc_flag + motor_flag) == 1);
assert(max(choice_flag + outcome_flag) == 1);


out.power_flags = {low_power_flag, high_power_flag};
out.area_flags = {frontal_flag, motor_flag, visual_flag, rsc_flag};
out.period_flags = {choice_flag, outcome_flag};


end