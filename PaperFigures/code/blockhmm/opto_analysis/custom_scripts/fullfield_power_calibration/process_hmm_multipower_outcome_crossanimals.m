%% plot hmm fractions and power dependence
epoch = 'outcome';
infof27 = get_hmm_info_animal('f27', epoch);
infof29 = get_hmm_info_animal('f29', epoch);
infof32 = get_hmm_info_animal('f32', epoch);

i = 6;
meanf27 = infof27.opto_frac_all(i,:) - infof27.opto_frac_all(i,1);

meanfracs = [meanf27(:, [1,3,4,5,6]);
    infof29.opto_frac_all(i,:) - infof29.opto_frac_all(i,1);
    infof32.opto_frac_all(i,:) - infof32.opto_frac_all(i,1)] * 100;


figure;
hold on
cols = brewermap(6, 'Set1');
cmap = cols([3,1,5,6,4,3],:);
cmap(6,:) = [0.5, 0.5, 0.5];

for i = 6
    errorbar(infof27.power_lst + i * 0.8, (infof27.opto_frac_all(i,:) - infof27.opto_frac_all(i,1)) * 100, infof27.opto_err_all(i,:) * 100 , 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 2)
    errorbar(infof29.power_lst + i * 0.8, meanfracs(2,:), infof29.opto_err_all(i,:) * 100, 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 2)
    errorbar(infof32.power_lst + i * 0.8, meanfracs(3,:), infof32.opto_err_all(i,:) * 100, 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 2)

    plot(infof32.power_lst + i * 0.8, mean(meanfracs, 1), 'o-', ...
        'Color', cmap(1,:), 'MarkerFaceColor', cmap(1,:), 'MarkerEdgeColor', cmap(1,:), ...
        'MarkerSize', 10, 'LineWidth', 2)
end

ylim([-50, 50])
hline(0, 'k--')
mymakeaxis('x_label', 'Power (mW)', 'y_label', 'Change in % of inference-based blocks', 'font_size', 25,...
    'xticks', 0:50:250, 'yticks', -50:50:50)


%% choice plot
epoch = 'choice';
infof27 = get_hmm_info_animal('f27', epoch);
infof29 = get_hmm_info_animal('f29', epoch);
% infof32 = get_hmm_info_animal('f32', epoch);

i = 6;
meanfracs = [infof27.opto_frac_all(i,:) - infof27.opto_frac_all(i,1);
    infof29.opto_frac_all(i,:) - infof29.opto_frac_all(i,1)] * 100;


figure;
hold on
cols = brewermap(6, 'Set1');
cmap = cols([3,1,5,6,4,3],:);
cmap(6,:) = [0.5, 0.5, 0.5];

for i = 6
    errorbar(infof27.power_lst + i * 0.8, meanfracs(1,:), infof27.opto_err_all(i,:) * 100 , 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 2)
    errorbar(infof29.power_lst + i * 0.8, meanfracs(2,:), infof29.opto_err_all(i,:) * 100, 'o-',...
        'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 2)

    plot(infof27.power_lst + i * 0.8, mean(meanfracs, 1), 'o-', ...
        'Color', cmap(1,:), 'MarkerFaceColor', cmap(1,:), 'MarkerEdgeColor', cmap(1,:), ...
        'MarkerSize', 10, 'LineWidth', 2)
end

ylim([-50, 50])
hline(0, 'k--')
mymakeaxis('x_label', 'Power (mW)', 'y_label', 'Change in % of inference-based blocks', 'font_size', 25,...
    'xticks', 0:50:250, 'yticks', -50:50:50)





function info = get_hmm_info_animal(animalname, epoch)
% Get the opto summary performance for given animal

if strcmp(epoch, 'outcome')
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
elseif strcmp(epoch, 'choice')
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

end

load(fullfile('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/opto_analysis/optodata/', expfitdate, 'opto_hmm_info.mat'));


info.opto_frac_all = nan(6, numel(power_lst));
info.noopto_frac_all = nan(6, numel(power_lst));
info.opto_err_all = nan(6, numel(power_lst));
info.noopto_err_all = nan(6, numel(power_lst));
info.power_lst = power_lst;


for i = 1:numel(power_lst)
    disp(power_lst(i))
    out = get_hmm_opto_frac(animalinfo, animalname, power_lst(i), tags{i}, epoch);
    assert(sum(out.sessid_lst) > 0)

    info.opto_frac_all(:,i) = out.opto_frac;
    info.opto_err_all(:,i) = out.opto_err;
    info.noopto_frac_all(:,i) = out.noopto_frac;
    info.noopto_err_all(:,i) = out.noopto_err;
    
end



end







function out = get_hmm_opto_frac(animalinfo, animal, power, tag, epoch)
% animalinfo: animalinfo object as output by the analysis script
% animalID: string, ID of animal to investigate
% power: power as an int
% tag: type of inactivation: 'OF': outcome with fake opto v3; 'O': outcome
% using previous version of opto inactivation protocol

animal_cands = strcmp({animalinfo.animal}, animal);
assert(sum(animal_cands) == 1)
animalID = find(animal_cands);

if strcmpi(epoch, 'outcome')
    if power == 240
        % for 240, include both 240 and 230
        sessstring = {sprintf('fullfield%d%s', power, lower(tag)), ...
            sprintf('fullfield%d%s', power - 10, lower(tag))};
    else
        sessstring = {sprintf('fullfield%d%s', power, lower(tag))};
    end
elseif strcmpi(epoch, 'choice')
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
else
    error('invalid epoch');
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


