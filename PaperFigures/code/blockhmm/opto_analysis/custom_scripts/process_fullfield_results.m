cd('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/custom_scripts')
% File for processing the animal info
% that was parsed by `hmmblockstate_visualize_022422_tsne_Nmodes.m`
% and `hmmstate_decode_single_animals_022422_tsne_Nmodes.m`
global NSTATES
NSTATES = 6;
expfitdate = '050422';

% Load the data
load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));



%% Extract fullfield/outcome sessions with high power based on criteria
animalID = 2; %f27 data
% sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O'); %179:181;
sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield120O');

zclasses = animalinfo(animalID).zclassified(sessid_lst);
opto = animalinfo(animalID).opto(sessid_lst);
obs = animalinfo(animalID).obs_splits(sessid_lst);
animal = animalinfo(animalID).animal;

count_opto = zeros(1, 6);
count_no_opto = zeros(1, 6);

% count number of blocks in each state
for id = 1:numel(opto)
    opto_sess = opto{id};
    zclass_sess = zclasses{id};

    for z = 1:6
        zcount_opto = sum(zclass_sess(opto_sess == 1) == z);
        zcount_no_opto = sum(zclass_sess(opto_sess == 0) == z);

        count_opto(z) = count_opto(z) + zcount_opto;
        count_no_opto(z) = count_no_opto(z) + zcount_no_opto;

    end
end

count_opto_f27 = count_opto;
count_no_opto_f27 = count_no_opto;
opto_frac_f27 = count_opto / sum(count_opto);
noopto_frac_f27 = count_no_opto / sum(count_no_opto);
opto_err_f27 = sqrt(opto_frac_f27 .* (1 - opto_frac_f27) ./ sum(count_opto));
noopto_err_f27 = sqrt(noopto_frac_f27 .* (1 - noopto_frac_f27) ./ sum(count_no_opto));


% 
% figure;
% plot(count_opto ./ sum(count_opto))
% hold on
% plot(count_no_opto ./ sum(count_no_opto))


%%
animalID = 3; %f29 data
% sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O');
sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield120O');


zclasses = animalinfo(animalID).zclassified(sessid_lst);
opto = animalinfo(animalID).opto(sessid_lst);
obs = animalinfo(animalID).obs_splits(sessid_lst);
animal = animalinfo(animalID).animal;

count_opto = zeros(1, 6);
count_no_opto = zeros(1, 6);

% count number of blocks in each state
for id = 1:numel(opto)
    opto_sess = opto{id};
    zclass_sess = zclasses{id};

    for z = 1:6
        zcount_opto = sum(zclass_sess(opto_sess == 1) == z);
        zcount_no_opto = sum(zclass_sess(opto_sess == 0) == z);

        count_opto(z) = count_opto(z) + zcount_opto;
        count_no_opto(z) = count_no_opto(z) + zcount_no_opto;

    end
end

count_opto_f29 = count_opto;
count_no_opto_f29 = count_no_opto;
opto_frac_f29 = count_opto / sum(count_opto);
noopto_frac_f29 = count_no_opto / sum(count_no_opto);
opto_err_f29 = sqrt(opto_frac_f29 .* (1 - opto_frac_f29) ./ sum(count_opto));
noopto_err_f29 = sqrt(noopto_frac_f29 .* (1 - noopto_frac_f29) ./ sum(count_no_opto));



%% Plot summary for f27/ f29
figure('Position', [440,485,782,313]);
msize = 10;
xvals_opto = (1:6) + 0.1;
xvals_no_opto = (1:6) - 0.1;
l3 = errorbar(xvals_opto, opto_frac_f27, opto_err_f27, 'bs',...
    'MarkerSize', msize, 'MarkerFaceColor', 'b');
hold on
l3b = errorbar(xvals_no_opto, noopto_frac_f27, noopto_err_f27, 'ks',...
    'MarkerSize', msize, 'MarkerFaceColor', 'k');

l4 = errorbar(xvals_opto, opto_frac_f29, opto_err_f29, 'bd',...
    'MarkerSize', msize, 'MarkerFaceColor', 'b');
hold on
l4b = errorbar(xvals_no_opto, noopto_frac_f29, noopto_err_f29, 'kd',...
    'MarkerSize', msize, 'MarkerFaceColor', 'k');

for i = 1:6
    plot([xvals_no_opto(i) xvals_opto(i)], [noopto_frac_f27(i), ...
            opto_frac_f27(i)], 'k');
    hold on
    plot([xvals_no_opto(i) xvals_opto(i)], [noopto_frac_f29(i), ...
            opto_frac_f29(i)], 'k');   
end
mymakeaxis('x_label', 'HMM States', 'y_label', 'Fraction', 'xticks', 1:6)
legend([l3, l3b, l4, l4b], {'f27 ON', 'f27 OFF', 'f29 ON', 'f29 OFF'}, ...
    'FontSize', 15)



%% Parse the transition functions
global NSTATES
NSTATES = 6;
expfitdate = '050422';
% animalID = 3; %f29 data
% sessid_lst = 176:178;
animalID = 2; %f27 data
% sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O');
sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield120O');


% Load the data
load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));
zclasses = animalinfo(animalID).zclassified(sessid_lst);
opto = animalinfo(animalID).opto(sessid_lst);
obs = animalinfo(animalID).obs_splits(sessid_lst);
animal = animalinfo(animalID).animal;

% Mean transition function
perf_opto = [];
perf_noopto = [];
for i = 1:numel(obs)
    opto_extract = obs{i}(opto{i} == 1, :);
    noopto_extract = obs{i}(opto{i} == 0, :);
    perf_opto = [perf_opto; opto_extract];
    perf_noopto = [perf_noopto; noopto_extract]; 
end

figure;
l1 = errorbar(1:size(perf_opto, 2), mean(perf_opto, 1), ...
    std(perf_opto, [], 1) / sqrt(size(perf_opto, 1)), 'b', 'LineWidth', 2);
hold on
l2 = errorbar(1:size(perf_noopto, 2), mean(perf_noopto, 1), ...
    std(perf_noopto, [], 1) / sqrt(size(perf_noopto, 1)), 'k', 'LineWidth', 2);
mymakeaxis('x_label', 'Trials in block', 'y_label', 'P(Correct)', ...
    'font_size', 25, 'xytitle', animal, 'xticks', 0:5:25)
legend([l1, l2], {'ON', 'OFF'}, 'FontSize', 15)



%%

power_criteria = {'low', 'high'};
area_criteria = {'frontal', 'visual', 'rsc', 'motor'};
period_criteria = {'outcome', 'choice'};


% Perform the parsing for all animals
modecounts = struct;
for animalID = 1:numel(animalinfo)
    fprintf('Processing animal %d of %d...\n', animalID, numel(animalinfo));
    out = helper.get_flags(animalinfo, animalID);
    
    
    zclasses = animalinfo(animalID).zclassified;
    opto = animalinfo(animalID).opto;
    obs = animalinfo(animalID).obs_splits;
    modecounts(animalID).animal = animalinfo(animalID).animal;
    
    count_opto_cell = cell(numel(power_criteria), numel(area_criteria),...
        numel(period_criteria));
    count_no_opto_cell = cell(numel(power_criteria), numel(area_criteria),...
        numel(period_criteria));
    transfunc_opto_cell = cell(numel(power_criteria), numel(area_criteria),...
        numel(period_criteria));
    transfunc_no_opto_cell = cell(numel(power_criteria), numel(area_criteria),...
        numel(period_criteria));
    
    for i = 1:numel(power_criteria)
        for j = 1:numel(area_criteria)
            for k = 1:numel(period_criteria)
                power_criterion = power_criteria{i};
                areas_criterion = area_criteria{j};
                period_criterion = period_criteria{k};
                
                [count_opto, count_no_opto, transfunc_opto, transfunc_no_opto] = ...
                    count_block_modes_from_criteria(zclasses, ...
                    opto, obs, out.power_flags, out.area_flags, ...
                    out.period_flags, power_criterion, areas_criterion, period_criterion);
                count_opto_cell{i,j,k} = count_opto;
                count_no_opto_cell{i,j,k} = count_no_opto;
                transfunc_opto_cell{i,j,k} = transfunc_opto;
                transfunc_no_opto_cell{i,j,k} = transfunc_no_opto;
%                 If we instead want to plot, uncomment this block of code
%                 plot_block_modes_from_criteria(zclasses, ...
%                     opto, out.power_flags, out.area_flags, ...
%                     out.period_flags, power_criterion, areas_criterion, period_criterion)
            end
        end
    end
    
    modecounts(animalID).count_opto = count_opto_cell;
    modecounts(animalID).count_no_opto = count_no_opto_cell;
    modecounts(animalID).power_criteria = power_criteria;
    modecounts(animalID).area_criteria = area_criteria;
    modecounts(animalID).period_criteria = period_criteria;
    modecounts(animalID).transfunc_opto = transfunc_opto_cell;
    modecounts(animalID).transfunc_no_opto = transfunc_no_opto_cell;
    modecounts(animalID).notes = 'count_opto(i,j,k) corresponds to power i, area j, period k';
end

% save
savefilename = fullfile('optodata/', expfitdate, 'modecounts.mat');


if exist(savefilename, 'file')
    response = questdlg(sprintf('File exists: %s, overwrite?', savefilename));
    switch response
        case 'Yes'
            save(savefilename, 'modecounts');
            fprintf('File overwritten\n');
        case 'No'
            fprintf('Skipped saving...\n');          
        case 'Cancel'
            error('Saving cancelled')
    end 
    
else
    save(savefilename, 'modecounts')
    fprintf('File saved!\n');
end


% 
% function out = get_flags(animalinfo, animalID)
% % animalinfo: struct saved by previous analysis
% % animalid: id from 1 - n of animal of interest
% % returns: out, structure containing all flags extracted from the animalinfo
% 
% % Parse information for the animal for each session
% opto_flag = cellfun(@(x) sum(x) > 0, animalinfo(animalID).opto);
% low_power_flag = cellfun(@(x) strcmp(x, 'Low'), animalinfo(animalID).power);
% high_power_flag = cellfun(@(x) strcmp(x, 'High'), animalinfo(animalID).power);
% frontal_flag = cellfun(@(x) contains(x, 'Frontal') & ~contains(x, 'Motor'), ...
%     animalinfo(animalID).areas);
% motor_flag = cellfun(@(x) contains(x, 'Motor') & ~contains(x, 'Visual'), ...
%     animalinfo(animalID).areas);
% rsc_flag = cellfun(@(x) contains(x, 'Rsc') & ~contains(x, 'Motor'), ...
%     animalinfo(animalID).areas);
% visual_flag = cellfun(@(x) contains(x, 'Visual') & ~contains(x, 'Motor'), ...
%     animalinfo(animalID).areas);
% choice_flag = cellfun(@(x) contains(x, 'Choice'), ...
%     animalinfo(animalID).period);
% outcome_flag = cellfun(@(x) contains(x, 'Outcome'), ...
%     animalinfo(animalID).period);
% 
% % make sure all flags are of the same size
% assert(max([numel(opto_flag), numel(low_power_flag), numel(frontal_flag),...
%     numel(motor_flag), numel(rsc_flag), numel(visual_flag)]) == numel(rsc_flag));
% assert(min([numel(opto_flag), numel(low_power_flag), numel(frontal_flag),...
%     numel(motor_flag), numel(rsc_flag), numel(visual_flag)]) == numel(rsc_flag));
% % assert(sum(low_power_flag + high_power_flag ~= opto_flag) == 0);
% % Now we will gather all blocks with the right opto information
% 
% out.power_flags = {low_power_flag, high_power_flag};
% out.area_flags = {frontal_flag, motor_flag, visual_flag, rsc_flag};
% out.period_flags = {choice_flag, outcome_flag};
% 
% 
% end



function plot_block_modes_from_criteria(zclasses, ...
    opto, power_flags, area_flags, ...
    period_flags, power_criterion, areas_criterion, period_criterion)
% zclasses: cell array of classes as extracted from animalinfo
% opto: cell array of opto locations as extracted from animalinfo
% power_flags: cell array containing high and low power flags
% area_flags: cell array containing flags for frontal, motor, visual and
% rsc
% period_flags: cell array containing choice/outcome flags
% power_criterion: 'low/high/none'
% areas_criterion: 'frontal/rsc/motor/visual/none'
% period_criterion: 'choice/outcome/none'
% if none: include all 
% produce two plots for counting the number and fraction of the HMM
% modes with and without laser, given the criterion
global NSTATES

[count_opto, count_no_opto] = count_block_modes_from_criteria(zclasses, ...
    opto, power_flags, area_flags, ...
    period_flags, power_criterion, areas_criterion, period_criterion);

% visualize
opto_fracs = count_opto / sum(count_opto);
no_opto_fracs = count_no_opto / sum(count_no_opto);
opto_err = sqrt(opto_fracs .* (1 - opto_fracs) / sum(count_opto));
no_opto_err = sqrt(no_opto_fracs .* (1 - no_opto_fracs) / sum(count_no_opto));

titlestring = sprintf('power = %s, area = %s, period = %s', power_criterion,...
    areas_criterion, period_criterion);

figure;
errorbar((1:NSTATES) - 0.1, opto_fracs, opto_err, 'bo', 'MarkerFaceColor', 'b')
hold on
errorbar((1:NSTATES) + 0.1, no_opto_fracs, no_opto_err, 'ro', 'MarkerFaceColor', 'r')
ylim([0, 1])
mymakeaxis('x_label', 'States', 'y_label', 'Fraction', 'xticks', 1:NSTATES,...
    'font_size', 25, 'xytitle', titlestring);


figure;
bar((1:NSTATES) - 0.1, count_opto, 'BarWidth', 0.1)
hold on
bar((1:NSTATES) + 0.1, count_no_opto, 'BarWidth', 0.1)
title(titlestring);
% ylim([0, 1])
% mymakeaxis('x_label', 'States', 'y_label', 'Fraction', 'xticks', 1:NSTATES,...
%     'font_size', 25);
end


function make_plot_states_and_opto(power_flags, area_flags, period_flags,...
    composition, animal_name)
%Power_flags: cell array containing high and low power flags
% area_flags: cell array containing flags for frontal, motor, visual and
% rsc
% period_flags: cell array containing choice/outcome flags
% composition: HMM mode composition
% will produce a standard plot of the HMM evolution together with
% information about which sessions are opto
global NSTATES

% Plot stacked bar plot and opto flags
figure('Position', [440,423,736,375], 'Name', animal_name);
h = bar(composition,'stacked', 'BarWidth', 1);
hold on
cols = paperaesthetics;
colors = cols.colorsLight;

for i = 1:NSTATES
    h(i).FaceColor = colors(i,:);
    h(i).ShowBaseLine = 'off';
    h(i).EdgeColor = 'none';
end



% Plot opto session with codes for brain area and power
sessions = 1:numel(power_flags{1});
levels_main = [0.5 0.8];
% sublevels = [-0.02 0.02];

for i = 1:numel(power_flags)
    for j = 1:numel(area_flags)
        for k = 1:numel(period_flags)
%             level = levels_main(i) + sublevels(k);
            level = levels_main(i);
            to_plot = sessions(power_flags{i} & area_flags{j} & ...
                period_flags{k});
            
            if numel(to_plot) == 0
                continue;
            end
            
            if k == 1
                marker = 'x';
            else
                marker = 'o';
            end
            
            switch j
                case 1
                    symbol = 'r';
                    l1 = plot(to_plot, level, [symbol marker], ...
                        'MarkerFaceColor', symbol);
                case 2
                    symbol = 'b';
                    l2 = plot(to_plot, level, [symbol marker], ...
                        'MarkerFaceColor', symbol);
                case 3
                    symbol = 'g';
                    l3 = plot(to_plot, level, [symbol marker], ...
                        'MarkerFaceColor', symbol);
                case 4
                    symbol = 'k';
                    l4 = plot(to_plot, level, [symbol marker], ...
                        'MarkerFaceColor', symbol);
                otherwise
                    error('Invalid area')
            end
            
            
           
        end
    end
end

ylim([0 1])
xlabel('Sessions')
ylabel('HMM state fraction')

% mymakeaxis('x_label', 'Sessions', 'y_label', 'Fraction', 'xytitle', ...
%     animal_name)
legend([l1(1), l2(1), l3(1), l4(1)], {'Motor', 'Vis.', 'Frontal', 'RSC'}, 'FontSize', 16);

set(gca, 'FontSize', 20);

end


function [count_opto, count_no_opto, transfunc_opto, transfunc_no_opto] = ...
    count_block_modes_from_criteria(zclasses, ...
    opto, obs, power_flags, area_flags, ...
    period_flags, power_criterion, areas_criterion, period_criterion)
% zclasses: cell array of classes as extracted from animalinfo
% opto: cell array of opto locations as extracted from animalinfo
% power_flags: cell array containing high and low power flags
% area_flags: cell array containing flags for frontal, motor, visual and
% rsc
% period_flags: cell array containing choice/outcome flags
% power_criterion: 'low/high/none'
% areas_criterion: 'frontal/rsc/motor/visual/none'
% period_criterion: 'choice/outcome/none'
% if none: include all 
% returns: count_opto: an 1x6 array with counts of opto blocks in each
% state
% count_no_opto: an 1x6 array with counts of non-opto blocks in each
% state
global NSTATES

[zclass_extracted, opto_extracted, obs_extracted, ~] = extract_block_modes_from_criteria(zclasses, ...
    opto, obs, power_flags, area_flags, ...
    period_flags, power_criterion, areas_criterion, period_criterion);
count_opto = nan(1, NSTATES);
count_no_opto = nan(1, NSTATES);

% count opto
zclass_opto = zclass_extracted(opto_extracted > 0);
obs_opto = obs_extracted(opto_extracted > 0, :);
for i = 1:NSTATES
    count_opto(i) = sum(zclass_opto == i);
%     transfunc_opto(i,:) = nanmean(obs_opto(zclass_opto == i, :), 1);
end
transfunc_opto = nanmean(obs_opto, 1);

% count non-opto
zclass_no_opto = zclass_extracted(opto_extracted == 0);
obs_no_opto = obs_extracted(opto_extracted == 0, :);

for i = 1:NSTATES
    count_no_opto(i) = sum(zclass_no_opto == i);
%     transfunc_no_opto(i,:) = nanmean(obs_no_opto(zclass_no_opto == i, :), 1);
end
transfunc_no_opto = nanmean(obs_no_opto, 1);



end




function [zclass_extracted, opto_extracted, obs_extracted, all_flags] = extract_block_modes_from_criteria(zclasses, ...
    opto, obs, power_flags, area_flags, ...
    period_flags, power_criterion, areas_criterion, period_criterion)
% zclasses: cell array of classes as extracted from animalinfo
% opto: cell array of opto locations as extracted from animalinfo
% power_flags: cell array containing high and low power flags
% area_flags: cell array containing flags for frontal, motor, visual and
% rsc
% period_flags: cell array containing choice/outcome flags
% power_criterion: 'low/high/none'
% areas_criterion: 'frontal/rsc/motor/visual/none'
% period_criterion: 'choice/outcome/none'
% if none: include all 
% Returns:
% zclass_extracted: 1 x nblocks array with each element representing a session,
% containing info about the z-classified modes in each block 
% opto_extracted: 1 x nblocks array with same size as zclass_extracted, 
% but giving information
% about the opto/non-opto identity of the block
% obs_extracted: cell array with each element representing a session, each
% element is nblocks x 15 array with outcomes in that block (starts at 0 at
% the beginning of each block and increases to 1).
% all_flags: array of size nsess x 1, indicating the sessions that are
% selected to extract

assert(numel(zclasses) == numel(opto));

switch power_criterion
    case 'low'
        power_flag = power_flags{1};
        
    case 'high'
        power_flag = power_flags{2};
        
    case 'none'
        power_flag = power_flags{1} + power_flags{2};
        
    otherwise
        error('invalid power')
end

switch period_criterion
    case 'choice'
        period_flag = period_flags{1};
        
    case 'outcome'
        period_flag = period_flags{2};
        
    case 'none'
        period_flag = period_flags{1} + period_flags{2};
        
    otherwise
        error('invalid power')
end


switch areas_criterion
    case 'frontal'
        area_flag = area_flags{1};
        
    case 'motor'
        area_flag = area_flags{2};
        
    case 'visual'
        area_flag = area_flags{3};
        
    case 'rsc'
        area_flag = area_flags{4};
        
    case 'none'
        area_flag = area_flags{1} + area_flags{2} + area_flags{3} + area_flags{4};
        
    otherwise
        error('invalid area')
end

assert(numel(area_flag) == numel(zclasses));
assert(numel(power_flag) == numel(zclasses));
assert(numel(period_flag) == numel(zclasses));


all_flags = area_flag & power_flag & period_flag;


zclass_extracted = cell2mat(zclasses(all_flags));
opto_extracted = cell2mat(opto(all_flags));
obs_extracted = cell2mat(obs(all_flags));

assert(size(obs_extracted, 1) == numel(opto_extracted));

if numel(opto_extracted) > 0
    assert(size(obs_extracted, 2) == 15);
end

end

