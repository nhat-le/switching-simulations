% For parsing the modecounts file generated by
% process_animal_info.m
expfitdate = '040122';

load(fullfile('optodata/', expfitdate, 'modecounts.mat'));

NSTATES = numel(modecounts(1).count_opto{1,1,1}); 

%% count the number of opto/no opto blocks for each condition per animal
nblocks_opto_all = cell(1, numel(modecounts));
nblocks_no_opto_all = cell(1, numel(modecounts));


for animalID = 1:numel(modecounts)
    count_opto = modecounts(animalID).count_opto;
    count_no_opto = modecounts(animalID).count_no_opto;
    
    nblocks_opto_all{animalID} = cellfun(@(x) sum(x), count_opto);
    nblocks_no_opto_all{animalID} = cellfun(@(x) sum(x), count_no_opto);
    
end

%%
period_criterion = 'choice';
nareas = 4;
msize = 4;
stateIDs = [2,3];
delta = 0.1; %how far apart opto/non-opto columns are
xvals_OFF = (1:nareas) - delta; %+ (animalID - 3) * delta / 3;
xvals_ON = (1:nareas) + delta; %+ (animalID - 3) * delta / 3;

figure('Position', [246,566,1114,232]);
[on1, off1] = plot_mode_fractions(modecounts, xvals_ON, xvals_OFF, [1], period_criterion, msize);

[on2, off2] = plot_mode_fractions(modecounts, xvals_ON + 5, xvals_OFF + 5, [2,3], period_criterion, msize);

[on3, off3] = plot_mode_fractions(modecounts, xvals_ON + 10, xvals_OFF + 10, [4], period_criterion, msize);


[on4, off4] = plot_mode_fractions(modecounts, xvals_ON + 15, xvals_OFF + 15, [5,6], period_criterion, msize);


region_labels = ioutils.capitalize(modecounts(1).area_criteria);
region_labels{3} = 'RSC';
region_labels = repmat(region_labels, [1, 4]);

ylim([0, 1])
mymakeaxis('x_label', 'Area', 'y_label', 'Fraction', 'font_size', 10, ...
    'xticks', [1:4, 6:9, 11:14, 16:19], 'xticklabels', region_labels, ...
    'yticks', [0, 0.5, 1])


% legend(lines, {'f26', 'f27', 'f29', 'f32'}, 'FontSize', 20);

%% Statistical tests
signrank(on4(:,3), off4(:,3))


%% plot of opto effect given specific conditions
% period_criterion = 'choice'; %oucome/choice/none
% power_criterion = 'low'; %low/high/none
% msize = 10;
% delta = 0.1; %how far apart opto/non-opto columns are
% 
% region_labels = ioutils.capitalize(modecounts(1).area_criteria);
% region_labels{3} = 'RSC';
% 
% 
% stateIDs = [2,3]; %1:NSTATES
% [idxperiod, idxpower] = helper.find_idx_set(modecounts, power_criterion, period_criterion);
% figure('Position', [440,478,682,320]);
% for animalID = 1:4 %numel(modecounts)
%     counts_opto = helper.block_counter(modecounts(animalID).count_opto, idxperiod, idxpower);
%     [nareas, nstates] = size(counts_opto);
% 
%     opto_states_frac = sum(counts_opto(:, stateIDs), 2) ./ sum(counts_opto, 2);
% 
% %     opto_states_frac = counts_opto ./ sum(counts_opto, 2);
%     opto_states_frac_error = sqrt(opto_states_frac .* (1 - opto_states_frac) ./ sum(counts_opto, 2));
% 
%     counts_no_opto = helper.block_counter(modecounts(animalID).count_no_opto, idxperiod, idxpower);
% %     no_opto_states_frac = counts_no_opto ./ sum(counts_no_opto, 2);
%     no_opto_states_frac = sum(counts_no_opto(:, stateIDs), 2) ./ sum(counts_no_opto, 2);
%     no_opto_states_frac_error = sqrt(no_opto_states_frac .* (1 - no_opto_states_frac) ./ sum(counts_no_opto, 2));
% 
% 
%     xvals_no_opto = (1:nareas) - delta + (animalID - 3) * delta / 3;
%     xvals_opto = (1:nareas) + delta + (animalID - 3) * delta / 3;
%     disp(xvals_no_opto)
% 
% 
%     % plot opto fractions
%     switch animalID        
%         case 1
%             l1 = errorbar(xvals_opto, opto_states_frac, opto_states_frac_error, 'bo',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_no_opto, no_opto_states_frac, no_opto_states_frac_error, 'ko',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');
% 
% 
%         case 2
%             l2 = errorbar(xvals_opto, opto_states_frac, opto_states_frac_error, 'bs',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_no_opto, no_opto_states_frac, no_opto_states_frac_error, 'ks',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');
% 
% 
% 
%         case 3
%             l3 = errorbar(xvals_opto, opto_states_frac, opto_states_frac_error, 'bd',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_no_opto, no_opto_states_frac, no_opto_states_frac_error, 'kd',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');
% 
% 
% 
%         case 4
%             l4 = errorbar(xvals_opto, opto_states_frac, opto_states_frac_error, 'b^',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_no_opto, no_opto_states_frac, no_opto_states_frac_error, 'k^',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');
% 
% 
% 
%     end     
%    % connect pairs of opto/no-opto
%     for areaID = 1:nareas
%         plot([xvals_no_opto(areaID), xvals_opto(areaID)], [no_opto_states_frac(areaID), ...
%             opto_states_frac(areaID)], 'k');   
%     end
% 


% end

% ylim([0, 1])
% mymakeaxis('x_label', 'Area', 'y_label', 'Fraction', 'font_size', 20, ...
%     'xticks', 1:4, 'xticklabels', region_labels)
% xticks(1:4)
% xticklabels(modecounts(1).area_criteria);
% 
% legend([l1, l2, l3, l4], {'f26', 'f27', 'f29', 'f32'}, 'FontSize', 20);



function [opto_frac_all, off_frac_all] = plot_mode_fractions(modecounts, xvals_ON, xvals_OFF, stateIDs, period_criterion, msize)
% xvals_ON: array of xvalues for ON condition
% xvals_OFF: array of xvalues for OFF condition
% stateIDs: the states to combine
% period_criterion: 'choice' or 'outcome'

% period_criterion = 'choice'; %oucome/choice/none
power_criterion = 'low'; %low/high/none
% msize = 10;
delta = 0.1; %how far apart opto/non-opto columns are

region_labels = ioutils.capitalize(modecounts(1).area_criteria);
region_labels{3} = 'RSC';


[idxperiod, idxpower] = helper.find_idx_set(modecounts, power_criterion, period_criterion);
opto_frac_all = [];
off_frac_all = [];


% figure('Position', [440,478,682,320]);
for animalID = 1:4 %numel(modecounts)
    counts_opto = helper.block_counter(modecounts(animalID).count_opto, idxperiod, idxpower);
    [nareas, ~] = size(counts_opto);

    opto_states_frac = sum(counts_opto(:, stateIDs), 2) ./ sum(counts_opto, 2);
    opto_frac_all(animalID, :) = opto_states_frac;

%     opto_states_frac = counts_opto ./ sum(counts_opto, 2);
    opto_states_frac_error = sqrt(opto_states_frac .* (1 - opto_states_frac) ./ sum(counts_opto, 2));

    counts_no_opto = helper.block_counter(modecounts(animalID).count_no_opto, idxperiod, idxpower);
%     no_opto_states_frac = counts_no_opto ./ sum(counts_no_opto, 2);
    no_opto_states_frac = sum(counts_no_opto(:, stateIDs), 2) ./ sum(counts_no_opto, 2);
    off_frac_all(animalID, :) = no_opto_states_frac;

    no_opto_states_frac_error = sqrt(no_opto_states_frac .* (1 - no_opto_states_frac) ./ sum(counts_no_opto, 2));

    xvals_OFF_plot = xvals_OFF; %+ (animalID - 3) * delta / 3;
    xvals_ON_plot = xvals_ON; %+ (animalID - 3) * delta / 3;
%     disp(xvals_OFF_plot);

%     xvals_OFF = (1:nareas) - delta + (animalID - 3) * delta / 3;
%     xvals_ON = (1:nareas) + delta + (animalID - 3) * delta / 3;


    % plot opto fractions
    switch animalID        
        case {1,2,3,4}
            % if we want error bars
%             l1 = errorbar(xvals_ON_plot, opto_states_frac, opto_states_frac_error, 'bo',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_OFF_plot, no_opto_states_frac, no_opto_states_frac_error, 'ko',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');

            l1 = plot(xvals_ON_plot, opto_states_frac, 'bo',...
                'MarkerSize', msize, 'MarkerFaceColor', 'b');
            hold on
            plot(xvals_OFF_plot, no_opto_states_frac, 'ko',...
                'MarkerSize', msize, 'MarkerFaceColor', 'k');


%         case 2
%             l2 = errorbar(xvals_ON_plot, opto_states_frac, opto_states_frac_error, 'bs',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_OFF_plot, no_opto_states_frac, no_opto_states_frac_error, 'ks',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');
% 
% 
% 
%         case 3
%             l3 = errorbar(xvals_ON_plot, opto_states_frac, opto_states_frac_error, 'bd',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_OFF_plot, no_opto_states_frac, no_opto_states_frac_error, 'kd',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');
% 
% 
% 
%         case 4
%             l4 = errorbar(xvals_ON_plot, opto_states_frac, opto_states_frac_error, 'b^',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'b');
%             hold on
%             errorbar(xvals_OFF_plot, no_opto_states_frac, no_opto_states_frac_error, 'k^',...
%                 'MarkerSize', msize, 'MarkerFaceColor', 'k');



    end     
   % connect pairs of opto/no-opto
    for areaID = 1:nareas
        plot([xvals_OFF_plot(areaID), xvals_ON(areaID)], [no_opto_states_frac(areaID), ...
            opto_states_frac(areaID)], 'k');   
    end



end

% lines = [l1, l2, l3, l4];

% ylim([0, 1])
% mymakeaxis('x_label', 'Area', 'y_label', 'Fraction', 'font_size', 20, ...
%     'xticks', 1:4, 'xticklabels', region_labels)
% xticks(1:4)
% xticklabels(modecounts(1).area_criteria);

% legend([l1, l2, l3, l4], {'f26', 'f27', 'f29', 'f32'}, 'FontSize', 20);


end
