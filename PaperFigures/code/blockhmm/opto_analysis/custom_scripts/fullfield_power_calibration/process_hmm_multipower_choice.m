cd('/Users/minhnhatle/Documents/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/custom_scripts')
% File for processing opto effect on HMM state composition,
% depending on the laser power

global NSTATES
NSTATES = 6;
animalname = 'f27';


switch animalname
    case 'f27'
        expfitdate = '052922';
        power_lst = [15, 30, 60, 120, 240]; 
        tags = {'C', 'C', 'C', 'C', 'C'};


    case 'f29'
        expfitdate = '051822';
        power_lst = [15, 30, 60, 120, 240]; 
        tags = {'C', 'C', 'C', 'C', 'C'};

    case 'f39'
        expfitdate = '080122';
        power_lst = [0, 15, 30, 60, 120, 240]; 
        tags = {'OF', 'OF', 'OF', 'OF', 'OF', 'OF'};

    otherwise
        error('Invalid animal')
end

load(fullfile('../optodata/', expfitdate, 'opto_hmm_info.mat'));


opto_frac_all = nan(6, numel(power_lst));
noopto_frac_all = nan(6, numel(power_lst));
opto_err_all = nan(6, numel(power_lst));
noopto_err_all = nan(6, numel(power_lst));


for i = 1:numel(power_lst)
    disp(power_lst(i))
    out = get_hmm_opto_frac(animalinfo, animalname, power_lst(i), tags{i});
    assert(sum(out.sessid_lst) > 0)

    opto_frac_all(:,i) = out.opto_frac;
    opto_err_all(:,i) = out.opto_err;
    noopto_frac_all(:,i) = out.noopto_frac;
    noopto_err_all(:,i) = out.noopto_err;
    
end

%% plot hmm fractions and power dependence
figure;
hold on
cols = brewermap(6, 'Set1');
cmap = cols([2,1,5,6,4,3],:);

for i = 1:6
    errorbar(power_lst + i * 0.8, opto_frac_all(i,:), opto_err_all(i,:), 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 10)
end

ylim([0, 1])
mymakeaxis('x_label', 'Power (mW)', 'y_label', 'Fraction', 'font_size', 25,...
    'xticks', 0:50:250, 'xytitle', animalname)


%% plot hmm fractions and power dependence (log scale)
% figure;
% hold on
% cols = brewermap(6, 'Set1');
% cmap = cols([2,1,5,6,4,3],:);
% 
% % i = 4;
% % plot(log(power_lst + 1) + i * 0.03, opto_frac_all(i,:), 'k', 'LineWidth', 2);
% for i = 1:6
%     errorbar(log10(power_lst + 1) + i * 0.01, opto_frac_all(i,:), opto_err_all(i,:), 'o-',...
%         'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
%         'MarkerSize', 10)
% end
% 
% 
% xlim([0, 3])
% ylim([0, 1])
% 
% 
% mymakeaxis('x_label', 'log Power (mW)', 'y_label', 'Fraction', 'font_size', 25,...
%     'xytitle', animalname)



function out = get_hmm_opto_frac(animalinfo, animal, power, tag)
% animalinfo: animalinfo object as output by the analysis script
% animalID: string, ID of animal to investigate
% power: power as an int
% tag: type of inactivation: 'OF': outcome with fake opto v3; 'O': outcome
% using previous version of opto inactivation protocol

animal_cands = strcmp({animalinfo.animal}, animal);
assert(sum(animal_cands) == 1)
animalID = find(animal_cands);

if strcmp(animal, 'f39')
    sessstring = {sprintf('fullfield%d%s', power, lower(tag))};
else
    if power == 240
        % for 240, include both 240 and 230
        sessstring = {sprintf('fullfield%s%d', lower(tag), power), ...
            sprintf('fullfield%d%s', power - 10, lower(tag))};
    else
        sessstring = {sprintf('fullfield%s%d', lower(tag), power)};
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


