paths = pathsetup('matchingsim');
% load(fullfile(paths.simdatapath, 'schematic/infbased_sim_combinations_122821.mat'));
load(fullfile(paths.simdatapath, 'deprecated/schematic/infbased_sim_combinations_122821c.mat'));


figure('Position', [432,554,703,244]);
hold on
colors = brewermap(9, 'Blues');
colors = colors([3, 5, 8], :);
type = 'epsilon'; %gamma or epsilon

coefs_arr = reshape(coefs_all, 10, 6, [], 3);

coefs_mean = squeeze(mean(coefs_arr, 1));
coefs_std = squeeze(std(coefs_arr, [], 1));


for id = 1:3
    subplot(1,3,id)
    
    if strcmp(type, 'gamma')
        idchoose = id;
    elseif strcmp(type, 'epsilon')
        idchoose = id + 3;
    else
        error('Invalid idchoose');
    end
    sub_choicearr = squeeze(choicelst_all(idchoose,2:101,:));
    imagesc(1-sub_choicearr, 'AlphaData', ~isnan(sub_choicearr));
    colormap redblue
    axis xy
    mymakeaxis('x_label', 'Trials from switch', 'y_label', 'Block #',...
        'xticks', 0:10:20)
    
%     subplot(4,3,id + 3)
%     errorbar(1:5, coefs_mean(id,:, 1), coefs_std(id,:,1), 'o-')
% %     ylim([-1 5])
%     hline(0)
%     
%     subplot(4,3,id + 6)
%     errorbar(1:5, coefs_mean(id,:, 2), coefs_std(id,:,2), 'o-')
% %     ylim([-1 15])
%     hline(0)
%     
%     subplot(4,3,id + 9)
%     errorbar(1:5, coefs_mean(id,:, 3), coefs_std(id,:,3), 'o-')
% %     ylim([-1 1])
%     hline(0)
    
end