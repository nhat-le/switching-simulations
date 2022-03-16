% Averaging results from two runs
folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/blockhmmfit/K_selection';
load(fullfile(folder, 'blockhmm_validation_021222.mat'));

ll_lst_all_copy = ll_lst_all;
load(fullfile(folder, 'blockhmm_validation_021322.mat'));
ll_lst_all = (ll_lst_all + ll_lst_all_copy) / 2;


figure('Position', [380,198,981,600]);
for i = 1:21
    subplot(3,7,i)
    plot(ll_lst_all(i,:));
%     title(animal_lst(i,:));
    ylim([0 4])
    xlim([0, 8])
    [~,idx] = max(ll_lst_all(i,:));
    vline(idx);
    
    mymakeaxis('x_label', 'K', 'y_label', 'c.v. LL', 'xytitle', animal_lst(i,:),...
        'xticks', [0, idx])
    
end