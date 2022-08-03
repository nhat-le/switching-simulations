%% Parse the transition functions
global NSTATES
NSTATES = 6;
clear sessid_lst
% expfitdate = '062322_f39';
expfitdate = '080122';
load(fullfile('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/optodata/', expfitdate, 'opto_hmm_info.mat'));

animalID = 5; 
% outcome inactivation
% sessid_lst1 = strcmpi(animalinfo(animalID).areas, 'Fullfield0OF');
% sessid_lst2 = strcmpi(animalinfo(animalID).areas, 'Fullfield15OF');
% sessid_lst3 = strcmpi(animalinfo(animalID).areas, 'Fullfield30OF');
% sessid_lst4 = strcmpi(animalinfo(animalID).areas, 'Fullfield60OF');
% sessid_lst5 = strcmpi(animalinfo(animalID).areas, 'Fullfield120OF');
% sessid_lst6 = strcmpi(animalinfo(animalID).areas, 'Fullfield230OF') | strcmpi(animalinfo(animalID).areas, 'Fullfield240OF');


% choice inactivation
sessid_lst1 = strcmpi(animalinfo(animalID).areas, 'Fullfield0CF');
sessid_lst2 = strcmpi(animalinfo(animalID).areas, 'Fullfield30CF');
sessid_lst3 = strcmpi(animalinfo(animalID).areas, 'Fullfield60CF');
sessid_lst4 = strcmpi(animalinfo(animalID).areas, 'Fullfield120CF');
sessid_lst5 = strcmpi(animalinfo(animalID).areas, 'Fullfield240CF');

sesslst_all = {sessid_lst1, sessid_lst2, sessid_lst3, sessid_lst4, sessid_lst5};
powers = [0 30 60 120 240];

% Load the data
animal = animalinfo(animalID).animal;

% Mean transition functions
perf_opto_all = {};
perf_noopto_all = {};
parr_all = {};
parr_means = [];
parr_stds = [];
Nbootstrap = 100;

for i = 1:numel(sesslst_all)
    if numel(sesslst_all{i} > 0)
        [perf_opto, perf_noopto] = get_perfs(animalinfo, animalID, sesslst_all{i});
        perf_opto_all{i} = perf_opto;
        perf_noopto_all{i} = perf_noopto;

        % Infer parameters
        parr = mathfuncs.fit_sigmoid_bootstrap(perf_opto, Nbootstrap);
        parr = parr(parr(:,3) > 0, :);

        parr_means(i,:) = nanmean(parr, 1);
        parr_stds(i,:) = nanstd(parr, [], 1);

        parr_all{i} = parr();


    else
        perf_opto_all{i} = nan;
        perf_noopto_all{i} = nan;
        parr_all{i} = nan;
        parr_means(i,:) = nan;
        parr_stds(i,:) = nan;
    end
end
%% collect and plot

figure;
figure;
subplot(131)
errorbar(powers, parr_means(:,1), parr_stds(:,1), 'o', 'Color', '#333191', ...
    'MarkerFaceColor', '#333191', 'MarkerEdgeColor', '#333191', 'LineWidth', 1);
xlim([0 250])
ylim([0, 3.5])
mymakeaxis('x_label', 'Power (mW)', 'xytitle', 'Offset', 'xticks', 0:100:200)

subplot(132)
errorbar(powers, parr_means(:,2), parr_stds(:,2), 'o', 'Color', '#333191', ...
    'MarkerFaceColor', '#333191', 'MarkerEdgeColor', '#333191', 'LineWidth', 1);
xlim([0 250])
ylim([0, 3.5])


mymakeaxis('x_label', 'Power (mW)', 'xytitle', 'Slope', 'xticks', 0:100:200)


subplot(133)
errorbar(powers, parr_means(:,3), parr_stds(:,3), 'o', 'Color', '#333191', ...
    'MarkerFaceColor', '#333191', 'MarkerEdgeColor', '#333191', 'LineWidth', 1);
xlim([0 250])
ylim([0, 0.4])


mymakeaxis('x_label', 'Power (mW)', 'xytitle', 'Lapse', 'xticks', 0:100:200)


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
    case '240'
        sessid_lst = strcmpi(animalinfo(animalID).areas, 'Fullfield240O') | strcmpi(animalinfo(animalID).areas, 'Fullfield230O');
    case '0'
        sessid_lst = strcmpi(animalinfo(animalID).areas, 'Fullfield0O');
    case '0F'
        sessid_lst = strcmpi(animalinfo(animalID).areas, 'Fullfield0OF');
    case '120'
        sessid_lst = strcmpi(animalinfo(animalID).areas, 'Fullfield120O');
    case '240F'
        sessid_lst = strcmpi(animalinfo(animalID).areas, 'Fullfield240OF')| strcmpi(animalinfo(animalID).areas, 'Fullfield230OF');
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


