%% Let's collect the information from all clusters

folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';
probs = [0 0.1 0.2 0.3];
meansAll = {};
stdsAll = {};
for i = 1:numel(probs)
    filename = get_filename(probs(i));
    load(fullfile(folder, filename));
    
    meansAll{i} = clusters.means;
    stdsAll{i} = clusters.stds;
end


%% Plotting
figure('Position', [626,401,721,353]);
hold on
paperaesthetics;
colors = brewermap(6, 'BuGn');
colors = colors(2:6,:);
xvals = linspace(-0.15, 0.15, 4);

mapvals = {[1 2 3 6], 1:6, [1 2 3 5 6], [1 2 3 5 6]};


plottype = {'Efficiency', 'Lapse', 'Slope', 'Offset'};
k = 4;

handles = [];
for i = 1:numel(probs)
    coltouse = colors(i,:);
    
    means = meansAll{i};
    if k == 4
        means = -means;
    end
    stds = stdsAll{i};
    h = errorbar(mapvals{i}' + xvals(i)', means(:,k), stds(:,k), 'o', 'MarkerFaceColor', coltouse,...
    'LineWidth', 2, 'MarkerEdgeColor', coltouse, 'Color', coltouse);
    handles(i) = h;
end

switch k
    case 1
        ylim([0.5 1])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:6, 'yticks', 0.5:0.1:1)

    case 2
        ylim([0 0.5])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:6, 'yticks', 0:0.1:0.5)

    case 3
        ylim([0 3])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:6, 'yticks', 0:1:3)
    
    case 4
        ylim([0 14])
        mymakeaxis('x_label', 'Class', 'y_label', plottype{k}, 'xticks', 1:6)
    
        
end
        
if k < 4
    l = legend(handles, {'1', '0.9', '0.8', '0.7'}, 'Position', [0.3766 0.7337 0.0749 0.1827]);
    l.Title.String = 'Probability';
    l.Color = 'none';
else
    l = legend(handles, {'1', '0.9', '0.8', '0.7'});
    l.Title.String = 'Probability';
    l.Color = 'none';
    
end


%% Decoding performance figure
arr = cellfun(@(x) find_perf(x, 2) * 100, counts_all);

%%
figure;
hold on
folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';
probs = [0 0.1 0.2 0.3];
mapvals = {[1 2 3 6], 1:6, [1 2 3 5 6], [1 2 3 5 6]};
colors = brewermap(6, 'BuGn');

handles = [];

for probi = 1:4
    filename = get_filename(probs(probi));
    load(fullfile(folder, filename));
    
    means = [];
    stds = [];
    coltouse = colors(probi + 2,:);
    for i = 1:numel(mapvals{probi})

        perf_clusti = cellfun(@(x) find_perf(x, i), counts_all);
        means(i) = mean(perf_clusti);
        stds(i) = std(perf_clusti);

    end


    h = errorbar(mapvals{probi}, means, stds, 'o', 'LineWidth', 1, ...
        'MarkerEdgeColor', coltouse, 'Color', coltouse);
    handles(probi) = h;
    plot(mapvals{probi}, means, 'Color', coltouse, 'LineWidth', 2);
%     plot(mapvals{probi}, means, 'LineWidth', 2);
end

mymakeaxis('x_label', 'Class', 'y_label', 'Decoding accuracy')
l = legend(handles, {'1', '0.9', '0.8', '0.7'});
l.FontSize = 12;
l.Title.String = 'Probability';
l.Title.FontSize = 12;
l.Color = 'none';


function res = find_perf(arr, i)
res = arr(i,i) / sum(arr(:,i));

end


function filename = get_filename(prob)
switch prob
    case 0
        filename = 'opts_prob0.0-2021-09-25 20.52.mat';
    case 0.1
        filename = 'opts_prob0.1-2021-09-25 21.44.mat';
    case 0.2
        filename = 'opts_prob0.2-2021-09-25 21.57.mat';
    case 0.3
        filename = 'opts_prob0.3-2021-09-25 22.29.mat';
end
end