cd('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/custom_scripts')
% File for processing the animal info
% that was parsed by `hmmblockstate_visualize_022422_tsne_Nmodes.m`
% and `hmmstate_decode_single_animals_022422_tsne_Nmodes.m`
global NSTATES
NSTATES = 6;
expfitdate = '062022';

% Load the data
load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));

%%
animals = {animalinfo.animal};
lst_outcome = {};
lst_choice = {};

for i = 1:numel(animals)
    fprintf('Animal: %s\n', animals{i})
    animal = animals{i};

    % outcome results
    regionstr = 'Fullfield240Of';
    period = 'Outcome8020';
    fprintf('** Gathering outcome sessions **\n')
    lst_outcome{i} = get_hmm_opto_frac(animalinfo, animal, regionstr, period);

    % choice results
    regionstr = 'Fullfield240Cf';
    period = 'Choice8020';
    fprintf('** Gathering choice sessions **\n')
    lst_choice{i} = get_hmm_opto_frac(animalinfo, animal, regionstr, period);
end


%% Plot the mode fractions
lst = lst_outcome;
figure('Position', [440,479,592,319]);
hold on

bstyles = {'bo', 'bs', 'bd', 'b^'};
kstyles = {'ko', 'ks', 'kd', 'k^'};


gap = 0.05;
handles = [];

for i = 1:numel(animals)
    msize = 10;
    xvals_opto = (1:6) + 0.1 + gap * i;
    xvals_no_opto = (1:6) - 0.1 + gap * i;
    
    handles(i) = errorbar(xvals_opto, lst{i}.opto_frac, lst{i}.opto_err, bstyles{i}, ...
        'MarkerSize', msize, 'MarkerFaceColor', 'b');
    errorbar(xvals_no_opto, lst{i}.noopto_frac, lst{i}.noopto_err, kstyles{i}, ...
        'MarkerSize', msize, 'MarkerFaceColor', 'k');
    
    % connect opto/no opto for same animal and modes
    for j = 1:6
        plot([xvals_no_opto(j), xvals_opto(j)], [lst{i}.noopto_frac(j), lst{i}.opto_frac(j)], 'k')

    end
    



end
mymakeaxis('x_label', 'HMM States', 'y_label', 'Fraction', 'xticks', 1:6)
legend(handles, animals, ...
    'FontSize', 15)

%% Plot the mode compositions as bar plots
% Get the compositions 
compositions_all = [];
lst = lst_outcome;
for i = 1:numel(lst)
    compositions_all = [compositions_all; lst{i}.noopto_frac; lst{i}.opto_frac];
end


figure('Position', [440,409,748,389]);
colors = brewermap(6, 'Set1');
colors = colors([2,1,5,6,4,3],:);
xvals = [1 2 4 5 7 8 10 11];
h = bar(xvals, compositions_all, 'stacked');

for i = 1:6
    h(i).FaceColor = colors(i,:);

end

mymakeaxis('x_label', '', 'y_label', 'HMM mode fraction', 'xytitle', 'Outcome 80-20 (fullfield)', ...
    'xticks', xvals, 'xticklabels', {'OFF', 'ON', 'OFF', 'ON', 'OFF', 'ON', 'OFF', 'ON'})


%% Plot the transition functions
lst = lst_choice;
pOpto_all = [];
pNoOpto_all = [];

savefolder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/figs/80-20_results/062022';

for i = 1:numel(animals)    
    % fit sigmoid to the trans functions
    transfunc_opto = mean(lst{i}.perf_opto, 1);
    transfunc_no_opto = mean(lst{i}.perf_noopto, 1);
    pOpto = mathfuncs.fit_sigmoid_asym(transfunc_opto);
    pNoOpto = mathfuncs.fit_sigmoid_asym(transfunc_no_opto);
    T = size(lst{i}.perf_opto, 2);

    pOpto_all(i, :) = pOpto;
    pNoOpto_all(i, :) = pNoOpto;

    figure;
    l1 = errorbar(1:T, transfunc_opto, ...
        std(lst{i}.perf_opto, [], 1) / sqrt(size(lst{i}.perf_opto, 1)), 'bo', 'LineWidth', 0.5,...
        'MarkerFaceColor', 'b');
    hold on
    l2 = errorbar(1:T, transfunc_no_opto, ...
        std(lst{i}.perf_noopto, [], 1) / sqrt(size(lst{i}.perf_noopto, 1)), 'ko', 'LineWidth', 0.5,...
        'MarkerFaceColor', 'k');

    titletext = animals{i};

    % plot curve fit
    plot(1:T, mathfuncs.sigmoid_asym(1:T, pOpto(1), pOpto(2), pOpto(3), pOpto(4)), 'b');
    plot(1:T, mathfuncs.sigmoid_asym(1:T, pNoOpto(1), pNoOpto(2), pNoOpto(3), pNoOpto(4)), 'k');

    ylim([0, 1])
    
    mymakeaxis('x_label', 'Trials in block', 'y_label', 'P (Correct)',...
        'xytitle', titletext, 'font_size', 20, 'xticks', 0:5:20);
    legend([l1, l2], {'Opto', 'No opto'}, 'FontSize', 20, 'Position', [0.698214285714286,0.404761911466205,0.203571428571429,0.125]);

    savefilename = sprintf('%s/transfunction-80-20-fullfield-240mW-%s-withfit-choice.pdf', ...
        savefolder, animals{i});
    saveas(gcf, savefilename);
    
end

%% Visualize the effect on the fit parameters
cols = paperaesthetics;
colid = 4; %1 = offset, 2 = slope, 3 = lapse
paramnames = {'Offset', 'Slope', 'LapseL', 'Lapse'};
figure('Position', [440,446,679,352]);

subplot(131)
hold on
plot_params(pOpto_all, pNoOpto_all, 1, paramnames)


subplot(132)
hold on
plot_params(pOpto_all, pNoOpto_all, 2, paramnames)


subplot(133)
hold on
plot_params(pOpto_all, pNoOpto_all, 4, paramnames)

%%

figure('Position', [440,446,679,352]);

subplot(131)
hold on

plot([-gap gap], max([pNoOpto_all(:,colid) pOpto_all(:,colid)], 0), 'k')
plot(ones(1, size(pNoOpto_all, 1)) * (-gap), max([pNoOpto_all(:,colid)], 0), 'ko', ...
    'MarkerFaceColor','k', 'MarkerEdgeColor', 'w', 'MarkerSize', 10)
plot(ones(1, size(pNoOpto_all, 1)) * (gap), max([pOpto_all(:,colid)], 0), 'bo', ...
    'MarkerFaceColor','b', 'MarkerEdgeColor', 'w', 'MarkerSize', 10)

mean1 = mean(pNoOpto_all(:,colid));
std1 = std(pNoOpto_all(:,colid));
mean2 = mean(pOpto_all(:,colid));
std2 = std(pOpto_all(:,colid));

plot([-gap gap], [mean1 mean2], 'Color', 'r', 'LineWidth', 2)


switch colid
    case {1,2}
        ylim([0 3])
    case {3,4}
        ylim([0 0.4])
end
xlim([-1, 1])

mymakeaxis('y_label', paramnames{colid},...
            'font_size', 25, ...
            'xticks', [-gap, gap], 'xticklabels', {'OFF', 'ON'});



function plot_params(pOpto_all, pNoOpto_all, colid, paramnames)
gap = 0.5;
plot([-gap gap], max([pNoOpto_all(:,colid) pOpto_all(:,colid)], 0), 'k')
plot(ones(1, size(pNoOpto_all, 1)) * (-gap), max([pNoOpto_all(:,colid)], 0), 'ko', ...
    'MarkerFaceColor','k', 'MarkerEdgeColor', 'w', 'MarkerSize', 10)
plot(ones(1, size(pNoOpto_all, 1)) * (gap), max([pOpto_all(:,colid)], 0), 'bo', ...
    'MarkerFaceColor','b', 'MarkerEdgeColor', 'w', 'MarkerSize', 10)

mean1 = mean(pNoOpto_all(:,colid));
std1 = std(pNoOpto_all(:,colid));
mean2 = mean(pOpto_all(:,colid));
std2 = std(pOpto_all(:,colid));

plot([-gap gap], [mean1 mean2], 'Color', 'r', 'LineWidth', 2)


switch colid
    case {1,2}
        ylim([0 3])
    case {3,4}
        ylim([0 0.6])
end
xlim([-1, 1])

mymakeaxis('y_label', paramnames{colid},...
            'font_size', 25, ...
            'xticks', [-gap, gap], 'xticklabels', {'OFF', 'ON'}, ...
            'customOffset', 0.5);
end


function out = get_hmm_opto_frac(animalinfo, animal, regionstr, period)
% animal: animal name
% sessid_lst: list of session ids to examine
% animalID: string, ID of animal to investigate
% returns: out, a structure with the counts of the z-classes and fractions
% together with the corresponding error bars

animal_cands = strcmp({animalinfo.animal}, animal);
assert(sum(animal_cands) == 1)
animalID = find(animal_cands);

sessid_lst = strcmp(animalinfo(animalID).areas, regionstr) & ...
    strcmp(animalinfo(animalID).period, period);
fprintf('Number of sessions = %d\n', sum(sessid_lst));


zclasses = animalinfo(animalID).zclassified(sessid_lst);

opto = animalinfo(animalID).opto(sessid_lst);
obs = animalinfo(animalID).obs_splits(sessid_lst);

out.animal = animalinfo(animalID).animal;

count_opto = zeros(1, 6);
count_no_opto = zeros(1, 6);

% count number of blocks in each state
for id = 1:numel(opto)
    opto_sess = opto{id};
    zclass_sess = zclasses{id};

    % TODO IMPORTANT: TO INVESTIGATE MODES 5 AND 6 FOR F32, TEMPORARY HACK
    if strcmp(animal, 'f32')
        zclass_sess(zclass_sess == 5) = 6;
    end

    for z = 1:6
        zcount_opto = sum(zclass_sess(opto_sess == 1) == z);
        zcount_no_opto = sum(zclass_sess(opto_sess == 0) == z);

        count_opto(z) = count_opto(z) + zcount_opto;
        count_no_opto(z) = count_no_opto(z) + zcount_no_opto;

    end
end

% Mean transition function
perf_opto = [];
perf_noopto = [];
for i = 1:numel(obs)
    opto_extract = obs{i}(opto{i} == 1, :);
    noopto_extract = obs{i}(opto{i} == 0, :);
    perf_opto = [perf_opto; opto_extract];
    perf_noopto = [perf_noopto; noopto_extract]; 
end


out.perf_opto = perf_opto;
out.perf_noopto = perf_noopto;
out.sessid_lst = sessid_lst;
out.count_opto = count_opto;
out.count_no_opto = count_no_opto;
out.opto_frac = count_opto / sum(count_opto);
out.noopto_frac = count_no_opto / sum(count_no_opto);
out.opto_err = sqrt(out.opto_frac .* (1 - out.opto_frac) ./ sum(count_opto));
out.noopto_err = sqrt(out.noopto_frac .* (1 - out.noopto_frac) ./ sum(count_no_opto));

end



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

