cd('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/custom_scripts')
% File for processing the animal info
% that was parsed by `hmmblockstate_visualize_022422_tsne_Nmodes.m`
% and `hmmstate_decode_single_animals_022422_tsne_Nmodes.m`
global NSTATES
NSTATES = 6;
expfitdate = '050822';

% Load the data
load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));



%% Extract fullfield/outcome sessions with high power based on criteria
% animalID = 2; %f27 data
% % sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O'); %179:181;
% sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield0O');
% 
% zclasses = animalinfo(animalID).zclassified(sessid_lst);
% opto = animalinfo(animalID).opto(sessid_lst);
% obs = animalinfo(animalID).obs_splits(sessid_lst);
% animal = animalinfo(animalID).animal;
% 
% count_opto = zeros(1, 6);
% count_no_opto = zeros(1, 6);
% 
% % count number of blocks in each state
% for id = 1:numel(opto)
%     opto_sess = opto{id};
%     zclass_sess = zclasses{id};
% 
%     for z = 1:6
%         zcount_opto = sum(zclass_sess(opto_sess == 1) == z);
%         zcount_no_opto = sum(zclass_sess(opto_sess == 0) == z);
% 
%         count_opto(z) = count_opto(z) + zcount_opto;
%         count_no_opto(z) = count_no_opto(z) + zcount_no_opto;
% 
%     end
% end
% 
% count_opto_f27 = count_opto;
% count_no_opto_f27 = count_no_opto;
% opto_frac_f27 = count_opto / sum(count_opto);
% noopto_frac_f27 = count_no_opto / sum(count_no_opto);
% opto_err_f27 = sqrt(opto_frac_f27 .* (1 - opto_frac_f27) ./ sum(count_opto));
% noopto_err_f27 = sqrt(noopto_frac_f27 .* (1 - noopto_frac_f27) ./ sum(count_no_opto));

%%
modetype = 0;
out_f27 = get_hmm_opto_frac(animalinfo, 2, modetype);
out_f29 = get_hmm_opto_frac(animalinfo, 3, modetype);

% 
% figure;
% plot(count_opto ./ sum(count_opto))
% hold on
% plot(count_no_opto ./ sum(count_no_opto))


%%
animalID = 3; %f29 data
% sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O');
sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield0O');


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
l3 = errorbar(xvals_opto, out_f27.opto_frac, out_f27.opto_err, 'bs',...
    'MarkerSize', msize, 'MarkerFaceColor', 'b');
hold on
l3b = errorbar(xvals_no_opto, out_f27.noopto_frac, out_f27.noopto_err, 'ks',...
    'MarkerSize', msize, 'MarkerFaceColor', 'k');

l4 = errorbar(xvals_opto, out_f29.opto_frac, out_f29.opto_err, 'bd',...
    'MarkerSize', msize, 'MarkerFaceColor', 'b');
hold on
l4b = errorbar(xvals_no_opto, out_f29.noopto_frac, out_f29.noopto_err, 'kd',...
    'MarkerSize', msize, 'MarkerFaceColor', 'k');

for i = 1:6
    plot([xvals_no_opto(i) xvals_opto(i)], [out_f27.noopto_frac(i), ...
            out_f27.opto_frac(i)], 'k');
    hold on
    plot([xvals_no_opto(i) xvals_opto(i)], [out_f29.noopto_frac(i), ...
            out_f29.opto_frac(i)], 'k');   
end
mymakeaxis('x_label', 'HMM States', 'y_label', 'Fraction', 'xticks', 1:6)
legend([l3, l3b, l4, l4b], {'f27 ON', 'f27 OFF', 'f29 ON', 'f29 OFF'}, ...
    'FontSize', 15)



%% Parse the transition functions
global NSTATES
NSTATES = 6;
clear sessid_lst
% expfitdate = '050422';
load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));

% animalID = 3; %f29 data
% sessid_lst = 176:178;
animalID = 2; %f27 data
sessid_lst1 = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O');
sessid_lst2 = strcmp(animalinfo(animalID).areas, 'Fullfield0O');

% Load the data
animal = animalinfo(animalID).animal;

% Mean transition functions
[perf_opto1, perf_noopto1] = get_perfs(animalinfo, animalID, sessid_lst1);
[perf_opto2, perf_noopto2] = get_perfs(animalinfo, animalID, sessid_lst2);


figure('Position', [24,162,554,407]);
l1 = errorbar(1:size(perf_opto1, 2), mean(perf_opto1, 1), ...
    std(perf_opto1, [], 1) / sqrt(size(perf_opto1, 1)), 'b', 'LineWidth', 2);
hold on
l2 = errorbar(1:size(perf_noopto1, 2), mean(perf_noopto1, 1), ...
    std(perf_noopto1, [], 1) / sqrt(size(perf_noopto1, 1)), 'k', 'LineWidth', 2);
l3 = errorbar(1:size(perf_opto2, 2), mean(perf_opto2, 1), ...
    std(perf_opto2, [], 1) / sqrt(size(perf_opto2, 1)), 'b--', 'LineWidth', 2);
l4 = errorbar(1:size(perf_noopto2, 2), mean(perf_noopto2, 1), ...
    std(perf_noopto2, [], 1) / sqrt(size(perf_noopto2, 1)), 'k--', 'LineWidth', 2);

mymakeaxis('x_label', 'Trials in block', 'y_label', 'P(Correct)', ...
    'font_size', 25, 'xytitle', animal, 'xticks', 0:5:25)
legend([l1, l2, l3, l4], {'ON 240', 'OFF 240', 'ON 0', 'OFF 0'}, 'FontSize', 15)



function [perf_opto, perf_noopto] = get_perfs(animalinfo, animalID, sessid_lst)
% animalinfo: animal info struct as saved by script
% animalID: index of animal to investigate
% sessid_lst: list of sessions to average over
% DEPRECATED: 
% zclasses: 1 x Nsess cell array, the z-state of HMM blocks in each session
% opto: 1 x Nsess cell array, opto or no-opto identity of each block in the
% session
% obs: Nsess x 1 cell array, each being Nblocks x T array: signed performance
% in each block of trials
opto = animalinfo(animalID).opto(sessid_lst);
obs = animalinfo(animalID).obs_splits(sessid_lst);

% Mean transition function
perf_opto = [];
perf_noopto = [];
for i = 1:numel(obs)
    opto_extract = obs{i}(opto{i} == 1, :);
    noopto_extract = obs{i}(opto{i} == 0, :);
    perf_opto = [perf_opto; opto_extract];
    perf_noopto = [perf_noopto; noopto_extract]; 
end


end


function out = get_hmm_opto_frac(animalinfo, animalID, modetype)
switch modetype
    case 240
        sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield240O') | strcmp(animalinfo(animalID).areas, 'Fullfield230O');
    case 0
        sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield0O');
    case 120
        sessid_lst = strcmp(animalinfo(animalID).areas, 'Fullfield0O');
    otherwise
        error('Invalid power')
end

zclasses = animalinfo(animalID).zclassified(sessid_lst);
opto = animalinfo(animalID).opto(sessid_lst);
out.animal = animalinfo(animalID).animal;

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

out.count_opto = count_opto;
out.count_no_opto = count_no_opto;
out.opto_frac = count_opto / sum(count_opto);
out.noopto_frac = count_no_opto / sum(count_no_opto);
out.opto_err = sqrt(out.opto_frac .* (1 - out.opto_frac) ./ sum(count_opto));
out.noopto_err = sqrt(out.noopto_frac .* (1 - out.noopto_frac) ./ sum(count_no_opto));

end


