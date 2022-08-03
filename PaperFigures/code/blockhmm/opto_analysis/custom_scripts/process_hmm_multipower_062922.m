cd('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/custom_scripts')
% File for processing opto effect on HMM state composition,
% depending on the laser power

global NSTATES
NSTATES = 6;

result_f27 = parse_opto_fracs('f27');
result_f29 = parse_opto_fracs('f29');
result_f32 = parse_opto_fracs('f32');

result_f27C = parse_opto_fracs_choice('f27');
result_f29C = parse_opto_fracs_choice('f29');


%% plot hmm fractions and power dependence
figure;
result = result_f32;
hold on
cols = brewermap(6, 'Set1');
cmap = cols([2,1,5,1,4,3],:);

for i = 4
    errorbar(result.power_lst + i * 0.8, result.opto_frac_all(i,:), result.opto_err_all(i,:), 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 10)

    errorbar(result.power_lst + i * 0.8, result.noopto_frac_all(i,:), result.noopto_err_all(i,:), 'o--',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 10)
end

ylim([0, 1])
mymakeaxis('x_label', 'Power (mW)', 'y_label', 'Fraction', 'font_size', 25,...
    'xticks', 0:50:250, 'xytitle', result.animal)


%% Opto effect single power

power = 120; %mw
results = {result_f27, result_f29, result_f32};
results_choice = {result_f27C, result_f29C};

modeIDs = [5, 6];
xloc = 1;
power = 120;
modeIDs = [6];

figure('Position', [440,478,892,320]);
plot_paired(1, results, power, [1]);
plot_paired(4, results, power, [2,3]);
plot_paired(7, results, power, [4]);
plot_paired(10, results, power, [5,6]);

% plot choice
plot_paired(16, results_choice, power, [1]);
plot_paired(19, results_choice, power, [2,3]);
plot_paired(22, results_choice, power, [4]);
plot_paired(25, results_choice, power, [5,6]);



labels = {'Mode 1', 'Mode 2+3', 'Mode 4', 'Mode 5+6', 'Mode 1', 'Mode 2+3', 'Mode 4', 'Mode 5+6'};

ylim([0, 1])
mymakeaxis('y_label', 'Fraction', 'font_size', 12, ...
    'xticks', [1.5 4.5 7.5 10.5 16.5, 19.5, 22.5, 25.5], 'xticklabels', labels, ...
    'yticks', [0, 0.5, 1])




%%
fracONs = [];
fracOFFs = [];
for i = 1:numel(results)
    result = results{i};
    idx = find(result.power_lst == power);

    opto_count = result.count_opto(:, idx);
    noopto_count = result.count_noopto(:, idx);

    fracON = sum(opto_count(modeIDs)) / sum(opto_count);
    fracOFF = sum(noopto_count(modeIDs)) / sum(noopto_count);

    
    fracONs(i) = fracON;
    fracOFFs(i) = fracOFF;

    plot([1,2], [fracOFF, fracON], 'k')
    hold on

end

plot(ones(1, numel(fracONs)) * 2, fracONs, 'bo', 'MarkerFaceColor', 'b')
plot(ones(1, numel(fracONs)), fracOFFs, 'ko', 'MarkerFaceColor', 'k')






function plot_paired(xloc, results, power, modeIDs)
% plot paired on/off fractions at the specified x-location
% inputs:
% - xloc: int, the location on the x axis to plot
% - results: result from the previous call
% - power: power in mW e.g. 60/120
% - modeIDs: list of modes to plot

fracONs = [];
fracOFFs = [];
for i = 1:numel(results)
    result = results{i};
    idx = find(result.power_lst == power);

    opto_count = result.count_opto(:, idx);
    noopto_count = result.count_noopto(:, idx);

    fracON = sum(opto_count(modeIDs)) / sum(opto_count);
    fracOFF = sum(noopto_count(modeIDs)) / sum(noopto_count);

    
    fracONs(i) = fracON;
    fracOFFs(i) = fracOFF;

    plot([xloc, xloc+1], [fracOFF, fracON], 'k')
    hold on

end

plot(ones(1, numel(fracONs)) * (xloc + 1), fracONs, 'bo', 'MarkerFaceColor', 'b')
plot(ones(1, numel(fracONs)) * xloc, fracOFFs, 'ko', 'MarkerFaceColor', 'k')


end




function result = parse_opto_fracs(animalname)
% Input: animalname: string
% Result: a struct containing the HMM fractions with or without opto

switch animalname
    case 'f27'
        expfitdate = '052922';
        power_lst = [0, 15, 30, 60, 120, 240]; 
        tags = {'OF', 'OF', 'OF', 'OF', 'O', 'OF'};


    case 'f29'
        expfitdate = '051822';
        power_lst = [0, 30, 60, 120, 240]; 
        tags = {'OF', 'OF', 'OF', 'O', 'OF'};

    case 'f32'
        expfitdate = '052922';
        power_lst = [0, 30, 60, 120, 240]; 
        tags = {'OF', 'OF', 'OF', 'OF', 'OF'};


    otherwise
        error('Invalid animal')
end

load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));

result.animal = animalname;
result.power_lst = power_lst;
result.opto_frac_all = nan(6, numel(power_lst));
result.noopto_frac_all = nan(6, numel(power_lst));
result.opto_err_all = nan(6, numel(power_lst));
result.noopto_err_all = nan(6, numel(power_lst));
result.count_opto = nan(6, numel(power_lst));
result.count_noopto = nan(6, numel(power_lst));




for i = 1:numel(power_lst)
    disp(power_lst(i))
    out = get_hmm_opto_frac(animalinfo, animalname, power_lst(i), tags{i});
    assert(sum(out.sessid_lst) > 0)

    result.opto_frac_all(:,i) = out.opto_frac;
    result.opto_err_all(:,i) = out.opto_err;
    result.noopto_frac_all(:,i) = out.noopto_frac;
    result.noopto_err_all(:,i) = out.noopto_err;
    result.count_opto(:,i) = out.count_opto;
    result.count_noopto(:,i) = out.count_no_opto;
    
end

end


function result = parse_opto_fracs_choice(animalname)
% Input: animalname: string
% Result: a struct containing the HMM fractions with or without opto

switch animalname
    case 'f27'
        expfitdate = '052922';
        power_lst = [0, 15, 30, 60, 120, 240]; 
        tags = {'OF', 'C', 'C', 'C', 'C', 'C'};


    case 'f29'
        expfitdate = '051822';
        power_lst = [0, 15, 30, 60, 120, 240]; 
        tags = {'OF', 'C', 'C', 'C', 'C', 'C'};

%     case 'f32'
%         expfitdate = '052922';
%         power_lst = [0, 30, 60, 120, 240]; 
%         tags = {'OF', 'OF', 'OF', 'OF', 'OF'};


    otherwise
        error('Invalid animal')
end

load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));

result.animal = animalname;
result.power_lst = power_lst;
result.opto_frac_all = nan(6, numel(power_lst));
result.noopto_frac_all = nan(6, numel(power_lst));
result.opto_err_all = nan(6, numel(power_lst));
result.noopto_err_all = nan(6, numel(power_lst));
result.count_opto = nan(6, numel(power_lst));
result.count_noopto = nan(6, numel(power_lst));




for i = 1:numel(power_lst)
    disp(power_lst(i))
    out = get_hmm_opto_frac(animalinfo, animalname, power_lst(i), tags{i});
    assert(sum(out.sessid_lst) > 0)

    result.opto_frac_all(:,i) = out.opto_frac;
    result.opto_err_all(:,i) = out.opto_err;
    result.noopto_frac_all(:,i) = out.noopto_frac;
    result.noopto_err_all(:,i) = out.noopto_err;
    result.count_opto(:,i) = out.count_opto;
    result.count_noopto(:,i) = out.count_no_opto;
    
end

end




function out = get_hmm_opto_frac(animalinfo, animal, power, tag)
% animalinfo: animalinfo object as output by the analysis script
% animalID: string, ID of animal to investigate
% power: power as an int
% tag: type of inactivation: 'OF': outcome with fake opto v3; 'O': outcome
% using previous version of opto inactivation protocol

animal_cands = strcmp({animalinfo.animal}, animal);
assert(sum(animal_cands) == 1)
animalID = find(animal_cands);

if power == 240
    if strcmpi(tag, 'C')
        sessstring = {sprintf('fullfieldc%d', power), ...
            sprintf('fullfieldc%d', power - 10)};
    else %outcome
    % for 240, include both 240 and 230
        sessstring = {sprintf('fullfield%d%s', power, lower(tag)), ...
            sprintf('fullfield%d%s', power - 10, lower(tag))};
    end
else
    if strcmpi(tag, 'C')
        sessstring = {sprintf('fullfieldc%d', power)};
    else %outcome
        sessstring = {sprintf('fullfield%d%s', power, lower(tag))};
    end
end

sessid_lst = contains(lower(animalinfo(animalID).areas), sessstring);

zclasses = animalinfo(animalID).zclassified(sessid_lst);




opto = animalinfo(animalID).opto(sessid_lst);
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

out.sessid_lst = sessid_lst;
out.count_opto = count_opto;
out.count_no_opto = count_no_opto;
out.opto_frac = count_opto / sum(count_opto);
out.noopto_frac = count_no_opto / sum(count_no_opto);
out.opto_err = sqrt(out.opto_frac .* (1 - out.opto_frac) ./ sum(count_opto));
out.noopto_err = sqrt(out.noopto_frac .* (1 - out.noopto_frac) ./ sum(count_no_opto));

end


