%% Parse the transition functions
global NSTATES
NSTATES = 6;
clear sessid_lst
expfitdate = '051822';
load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));

animalID = 2; %f27 data
% outcome inactivation
sessid_lst1 = strcmpi(animalinfo(animalID).areas, 'Fullfield0OF');
sessid_lst2 = strcmpi(animalinfo(animalID).areas, 'Fullfield30OF');
sessid_lst3 = strcmpi(animalinfo(animalID).areas, 'Fullfield60OF');
sessid_lst4 = strcmpi(animalinfo(animalID).areas, 'Fullfield120O');
sessid_lst5 = strcmpi(animalinfo(animalID).areas, 'Fullfield230OF') | strcmpi(animalinfo(animalID).areas, 'Fullfield240OF');


% choice inactivation
% sessid_lst1 = strcmpi(animalinfo(animalID).areas, 'FullfieldC60');
% sessid_lst2 = strcmpi(animalinfo(animalID).areas, 'FullfieldC120');
% sessid_lst3 = strcmpi(animalinfo(animalID).areas, 'FullfieldC240');
% sessid_lst4 = []; %strcmpi(animalinfo(animalID).areas, 'Fullfield230OF') | strcmpi(animalinfo(animalID).areas, 'Fullfield240OF');

sesslst_all = {sessid_lst1, sessid_lst2, sessid_lst3, sessid_lst4, sessid_lst5};
powers = [0 30 60 120 240];



%%


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

        parr_means(i,:) = mean(parr, 1);
        parr_stds(i,:) = std(parr, [], 1);

        parr_all{i} = parr();


    else
        perf_opto_all{i} = nan;
        perf_noopto_all{i} = nan;
        parr_all{i} = nan;
    end
end
%% collect and plot

figure;
subplot(131)
errorbar(powers, parr_means(:,1), parr_stds(:,1), 'o');
xlim([0 250])
mymakeaxis('x_label', 'Power (mW)', 'xytitle', 'Offset', 'xticks', 0:100:200)

subplot(132)
errorbar(powers, parr_means(:,2), parr_stds(:,2), 'o');
xlim([0 250])

mymakeaxis('x_label', 'Power (mW)', 'xytitle', 'Slope', 'xticks', 0:100:200)


subplot(133)
errorbar(powers, parr_means(:,3), parr_stds(:,3), 'o');
xlim([0 250])

mymakeaxis('x_label', 'Power (mW)', 'xytitle', 'Lapse', 'xticks', 0:100:200)



%%
pOpto1 = mathfuncs.fit_sigmoid(mean(perf_opto1));
pOpto2 = mathfuncs.fit_sigmoid(mean(perf_opto2));
pOpto3 = mathfuncs.fit_sigmoid(mean(perf_opto3));
pOpto4 = mathfuncs.fit_sigmoid(mean(perf_opto4));

pOff1 = mathfuncs.fit_sigmoid(mean(perf_noopto1));
pOff2 = mathfuncs.fit_sigmoid(mean(perf_noopto2));
pOff3 = mathfuncs.fit_sigmoid(mean(perf_noopto3));
pOff4 = mathfuncs.fit_sigmoid(mean(perf_noopto4));

powers = [0 60 120 240];

figure;
subplot(131)
plot(powers, [pOpto1(1) pOpto2(1) pOpto3(1) pOpto4(1)])
hold on
plot(powers, [pOff1(1) pOff2(1) pOff3(1) pOff4(1)])
mymakeaxis('x_label', 'Power', 'xytitle', 'Offset')


subplot(132)
plot(powers, [pOpto1(2) pOpto2(2) pOpto3(2) pOpto4(2)])
hold on
plot(powers, [pOff1(2) pOff2(2) pOff3(2) pOff4(2)])
mymakeaxis('x_label', 'Power', 'xytitle', 'Slope')

subplot(133)
plot(powers, [pOpto1(3) pOpto2(3) pOpto3(3) pOpto4(3)])
hold on
plot(powers, [pOff1(3) pOff2(3) pOff3(3) pOff4(3)])
mymakeaxis('x_label', 'Power', 'xytitle', 'Lapse')



% xvals = 1:15;
% plot(xvals, mathfuncs.sigmoid(xvals, pOpto(1), pOpto(2), pOpto(3)));
% hold on;
% plot(meanperf, 'o')







%%

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
% legend([l1, l2, l3, l4], {'ON 240v2', 'OFF 240v2', 'ON 240v3', 'OFF 240v3'}, 'FontSize', 15)
legend([l3, l4], {'ON 60mW', 'OFF 60mW'}, 'FontSize', 15)



function L = findLLH(params, y)
    % likelihood function of the observations y
    % given parameters params
    mu = params(1); % Center of transition
    sigma = params(2); % slope of transition
    eps = params(3);


    xvals = 1:numel(y);
    p = 1 ./ (1 + exp(-((xvals - mu) / sigma)));
    p = min(p, 0.9999);
    p = max(p, 0.0001);
    
    % Discard locations where y==0.5 (miss)
    p = p(y ~= 0.5);
    y = y(y ~= 0.5);

    L = - sum((1-y) .* log(1-p) + y .* log(p));

end



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


