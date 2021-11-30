% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/svm_Mdl.mat');
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params.mat');

% 9.29.21 expfit params updated
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_092921.mat');

% 9.30.21 prob expfit params
% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params_prob_093021.mat');

% load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_100821_withknnMdl.mat')
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_101421_withknnMdl.mat')


Mdl = Mdls1{1};
% Decode for all animal sessions
lapseFlat = reshape(explapses_all, [], 1);
effFlat = reshape(expeff_all, [], 1);
offsetFlat = reshape(expoffsets_all, [], 1);
slopesFlat = reshape(expslopes_all, [], 1);

% Note: no normalization since normalization already handled in svm Mdl
features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];

statesFlat = Mdl.predict(features_flat);
% eliminate nan's (nan's decode to 1)
statesFlat(isnan(lapseFlat) | isnan(effFlat) | isnan(offsetFlat) | isnan(slopesFlat)) = nan;


% Unflatten
states = reshape(statesFlat, size(expeff_all, 1), []);

%%
statecounts = [];
for i = 1:5
    statecounts(i,:) = sum(states == i, 1);
end

colors = brewermap(6, 'Set1');

statefrac = statecounts ./ sum(statecounts);
statecumulative = cumsum(statefrac);

figure('Position', [440,423,736,375]);
h = bar(statefrac','stacked');
orders = [2, 1, 5, 4, 3];
for i = 1:5
    h(i).FaceColor = colors(orders(i),:);
    h(i).ShowBaseLine = 'off';
end
xlim([0.5,16.5])
mymakeaxis('x_label', 'Training days', 'y_label', 'Fraction of animals', 'xticks', 1:16)
l = legend(h(1:3), {'Class 1', 'Class 2', 'Class 3'}, 'Position', [0.4,0.42,0.1,0.1], ...
    'FontSize', 14);

%%
paperaesthetics;
errorbar(1:size(states, 2), nanmean(states), nanstd(states, [], 1) / sqrt(size(states, 1)), 'o',...
    'MarkerFaceColor', bluecol, 'MarkerEdgeColor', bluecol, 'Color', bluecol,...
    'LineWidth', 2)
xlim([0, 16])
ylim([1 3])
mymakeaxis('x_label', 'Training days', 'y_label', 'Mean decoded state', 'yticks', 1:3)


%% Visualize animal behavior in mds space 
usepca = 0;
[res1, ~] = load_and_run(0, 1);
[res3, ~] = load_and_run(0.2);

% Split/combine so that we end up with 5 clusters for each space
% For prob = 1, will split cluster 4 into two, based on results of clusters
% for prob = 0.8
res1.idx(res1.idx == 4 & res3.idx == 5) = 5;
labels = res1.idx;


features_flatcopy = features_flat;
features_flatcopy(:,3) = -features_flat(:,3);
features_flatcopy(features_flatcopy < -20) = 3;

% features_combined = [features_flatcopy; res1.features];
features_norm = (features_flatcopy - nanmean(res1.features, 1)) ./ nanstd(res1.features, [], 1);
features_proj = features_norm * res1.V;
Yanimals = reshape(features_proj(:, 1:2), size(expeff_all, 1), [], 2);


%%
% focus on f12 for main figure
colors = brewermap(6, 'Pastel1');
colorsbold = brewermap(6, 'Set1');

colors = colors([2,1,5,4,3],:);
colorsbold = colorsbold([2,1,5,4,3],:);


figure;
hold on

% nearest-neighbor decoding
Y = res1.Y;
xmin = -6.5;
xmax = 4;
ymin = -4;
ymax = 4;

xcoords = linspace(xmin, xmax, 100);
ycoords = linspace(ymin, ymax, 100);
% [xx, yy] = meshgrid(xcoords, ycoords, ycoords);

domains = zeros(numel(xcoords), numel(ycoords));

% for each (x,y), find the Y entry that is closest
for i = 1:numel(xcoords)
    for j = 1:numel(ycoords)
        xval = xcoords(i);
        yval = ycoords(j);
        
        D = [xval yval] - res1.Y;
        
        dist = D(:,1).^2 + D(:,2).^2;
        
        domains(i,j) = Mdls1{1}.predict([xval; yval]);
        
%         id = argmin(dist);
        
%         domains(i,j) = res1.idx(id);
        
        
    end
end

imagesc(domains', 'XData', xcoords, 'YData', ycoords);
for i = 1:5
    plot(res1.Y(res1.idx == i, 1), res1.Y(res1.idx == i, 2), 'o', ...
      'MarkerFaceColor', colorsbold(i,:), 'MarkerEdgeColor', 'w', 'MarkerSize', 7);

end
i = 5;
plot(Yanimals(i,1,1), Yanimals(i,1,2), 'rx', 'LineWidth', 3, 'MarkerSize', 15)
plot(Yanimals(i,:, 1), Yanimals(i,:,2), 'k', 'LineWidth', 2)
colormap(colors)
% title(all_animals(i,:))

xlim([xmin, xmax])
ylim([ymin, ymax])
mymakeaxis('x_label', 'PC1', 'y_label', 'PC2')


%% plot all animals
figure;
hold on
for i = 1:15
    subplot(4,4,i)
    hold on
    
    imagesc(domains', 'XData', xcoords, 'YData', ycoords);
%     for j = 1:5
%         plot(res1.Y(res1.idx == j, 1), res1.Y(res1.idx == j, 2), 'o', ...
%           'MarkerFaceColor', colorsbold(j,:), 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
% 
%     end
    plot(Yanimals(i,1,1), Yanimals(i,1,2), 'rx', 'LineWidth', 3, 'MarkerSize', 15)
    plot(Yanimals(i,:, 1), Yanimals(i,:,2), 'k', 'LineWidth', 2)
    colormap(colors)
    

    xlim([xmin, xmax])
    ylim([ymin, ymax])
    
    
    title(all_animals(i,:))
    set(gca, 'FontSize', 10, 'FontName', 'helvetica');
end



