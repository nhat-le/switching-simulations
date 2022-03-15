% Code for fitting rho for all behavioral sessions
% Will use the fitted rho data in the directory
% processed_data/rho_data/{version}

paths = pathsetup('matchingsim');
version = '122221b';

filedir = fullfile(paths.datapath, 'rho_data', version, '*.mat');
files = dir(filedir);

rho_cell = {};

for i = 1:numel(files)
    load(fullfile(files(i).folder, files(i).name))
    rho_cell{i} = rhos;
end


rho_arr = pad_to_same_length(rho_cell);

rho_means = nanmean(rho_arr, 1);
rho_stds = nanstd(rho_arr, [], 1);
sesscounts = sum(~isnan(rho_arr), 1);
rho_stderr = rho_stds ./ sqrt(sesscounts);
display_range = 1:25;
figure;
errorbar(display_range, rho_means(display_range), rho_stderr(display_range), 'o',...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
hold on
hline(0, 'k--')
mymakeaxis('x_label', 'Session', 'y_label', '\rho',...
    'font_size', 30)


%% Statistical tests
% first five sessions > 0?
rho_first_five = reshape(rho_arr(:,1:5), 1, []);
rho_first_five = rho_first_five(~isnan(rho_first_five));
p1 = signrank(rho_first_five);
%p1 < 1e-5


% last five sessions > 0?
rho_last_five = reshape(rho_arr(:,21:25), 1, []);
rho_last_five = rho_last_five(~isnan(rho_last_five));
p2 = signrank(rho_last_five);
% p2 = 0.32


%% Individual animal plots
paths = pathsetup('matchingsim');
cols = paperaesthetics;

load(fullfile(paths.datapath, 'rho_data/122221b/f11_rhofit_122221b.mat'));
[rhof11, sef11] = find_rho_and_se(Ne_lst, Nr_lst);

load(fullfile(paths.datapath, 'rho_data/122221b/f16_rhofit_122221b.mat'));
[rhof16, sef16] = find_rho_and_se(Ne_lst, Nr_lst);


figure;
% l1 = plot(rho_arr(11,1:40)', 'b'); % f11
l1 = errorbar(1:numel(rhof11), rhof11, sef11, 'Color', cols.redcol);
hold on
% plot(rho_arr(13,1:40)', 'r') % f16
% l2 = plot(rho_arr(3,1:40)', 'r'); % e53
l2 = errorbar(1:numel(rhof16), rhof16, sef16, 'Color', cols.bluecol);
xlim([1, 25])


mymakeaxis('x_label', 'Session', 'y_label', '\rho',...
    'font_size', 30, 'yticks', -1:0.5:1, 'xticks', 0:5:25)

legend([l1, l2], {'f11', 'f16'}, 'FontSize', 20)


%% Plot rho evolution for all animals
paths = pathsetup('matchingsim');
version = '122221b';
filedir = fullfile(paths.datapath, 'rho_data', version, '*.mat');
files = dir(filedir);

for i = 1 :numel(files)
    load(fullfile(files(i).folder, files(i).name))
    parts = strsplit(files(i).name, '_');
    animal_name = parts{1};
    
    [rho, se] = find_rho_and_se(Ne_lst, Nr_lst);
    
    figure;
    errorbar(1:numel(rho), rho, se, 'Color', 'k');
    hold on
    hline(0, 'k--');
    xmax = min(numel(rho), 37);
    xlim([0, xmax]);
    ylim([-1, 1.2])


    mymakeaxis('x_label', 'Session', 'y_label', '\rho',...
        'font_size', 30, 'yticks', -1:0.5:1, 'xticks', 0:5:xmax,...
        'xytitle', animal_name)
    

    
    
    savefilename = fullfile(paths.figpath, sprintf('rho_expfit/%s_rhoexpfit.pdf', animal_name));
    if ~exist(savefilename, 'file')
        saveas(gcf, savefilename);
        fprintf('File saved!\n');
    end
end




function [rho, se] = find_rho_and_se(Ne_lst, Nr_lst)
% Given Ne_lst and Nr_lst: cells of Ne and Nr counts
% for individual session
% Returns: rho: array of estimated rho values
% se: standard error of the estimates
assert(numel(Ne_lst) == numel(Nr_lst))
rho = [];
for i = 1:numel(Ne_lst)
    Ne = double(Ne_lst{i}');
    Nr = double(Nr_lst{i}');
    
    if numel(Ne) == 0
        rho(i) = nan;
        se(i) = nan;
    else
    
        rho(i) = corr(Ne, Nr);
        se(i) = sqrt((1 - rho(i)^2) / (numel(Ne) - 2));
    end
end




end






