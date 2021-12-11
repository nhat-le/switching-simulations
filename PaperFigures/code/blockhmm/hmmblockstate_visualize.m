
rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/';
folders = dir(fullfile(rootdir, '*hmmblockfit_102121.mat'));

aggmeans = {};
aggparams = {};
for i = 1:numel(folders)
    load(fullfile(rootdir, folders(i).name));
    
    if i == numel(folders)
        obs(isnan(obs)) = 1;
    end
    
    allmeans = getmeans(obs, zstates);
    aggmeans{i} = allmeans;
    aggparams{i} = params;
     
    
end

aggmeansarr = cell2mat(aggmeans');
aggparamsarr = cell2mat(aggparams');


% try to refit the sigmoidal for all the curves

%%
ytarget = allmeans(1,:);
xvals = 0:14;
yguess = sigmoid([4, 2, 0.1], xvals);

pfit = fminsearch(@(p) sigmoidloss(p, xvals, ytarget), [4, 2, 0.1]);
ypred = sigmoid(pfit, xvals);
figure;
plot(xvals, ytarget);
hold on
plot(xvals, ypred);


%%
pfitall = [];
ypredall = [];
for i = 1:size(aggmeansarr, 1)
    ytarget = aggmeansarr(i,:);
    pfit = fminsearch(@(p) sigmoidloss(p, xvals, ytarget), [4, 2, 0.1]);
    pfitall(i,:) = pfit; 
    ypred = sigmoid(pfit, xvals);
    ypredall(i,:) = ypred;

end

%%
figure;
for i = 1:48
    subplot(12,4,i);
    plot(xvals, aggmeansarr(i,:));
    hold on
    plot(xvals, ypredall(i,:));
 
end


%%
titles = {'mu', 'sigma', 'lapse', 'unknown'};
figure;
for i = 1:4
    for j = 1:4
        subplot(4,4,4*(i-1) + j);
        plot(aggparamsarr(:,i), aggparamsarr(:,j), '.')
        xlim([0, 3])
        ylim([0,3])
        title(titles{i});
    end
end


%%
figure;
titles = {'mu', 'sigma', 'lapse'};

for i = 1:3
    for j = 1:3
        subplot(3,3,3*(i-1) + j);
        plot(pfitall(:,i), pfitall(:,j), '.')
        xlim([0, 3])
        ylim([0,3])
        title(titles{i});
    end
end

%% manual scoring
manualclass = [2, 3, 3, 1;
    1, 2, 3, 4;
    3, 4, 2, 4;
    1, 3, 3, 2;
    3, 3, 3, 3;
    4, 1, 3, 2;
    4, 2, 1, 1;
    3, 2, 3, 3;
    1, 3, 2, 4;
    3, 1, 2, 3;
    3, 3, 2, 3;
    2, 3, 3, 2;
    2, 4, 1, 3; %end1
    1, 2, 1, 3;
    3, 2, 2, 3;
    2, 1, 2, 4; %mode 4 unsure?
    2, 3, 3, 2;
    3, 2, 4, 1;
    3, 3, 3, 2;
    1, 3, 2, 3;
    3, 3, 2, 1;
    2, 3, 1, 4;
    4, 2, 3, 4];

manualflat = reshape(manualclass', 1, []);

classflat3 = (pfitall(:,1) < 2.5) & (pfitall(:,3) < 0.15); %inf based class
classflat2 = (pfitall(:,1) >= 2.5) & (pfitall(:,1) < 10) & (pfitall(:,3) < 0.2) & (pfitall(:,2) < 2);
% qlearning class
% classflat3 = (pfitall(:,1) >= 10);
classflat4 = (pfitall(:,1) < 2.5) & (pfitall(:,3) >= 0.15);
% reversal class

classflat1 = ~(classflat4 | classflat2 | classflat3);
%random class

classflat = classflat1 + classflat2 * 2 + classflat3 * 3 + classflat4 * 4;
classarr = reshape(classflat, 4, [])';

%% save
rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata';
% save(fullfile(rootdir, 'allanimals_hmmblockfitsummary_102321.mat'), 'aggparamsarr', 'classarr', 'pfitall',...
%     'ypredall', 'aggmeansarr', 'classflat', 'manualflat', 'manualclass', 'folders')


%%
pfitsall1 = pfitall(manualflat == 1,:);
pfitsall2 = pfitall(manualflat == 2,:);
pfitsall3 = pfitall(manualflat == 3,:);
pfitsall4 = pfitall(manualflat == 4,:);

figure;
subplot(131)
plot(pfitsall1(:,1), pfitsall1(:,2), 'x');
hold on
plot(pfitsall2(:,1), pfitsall2(:,2), 'x');
plot(pfitsall3(:,1), pfitsall3(:,2), 'x');
plot(pfitsall4(:,1), pfitsall4(:,2), 'x');
xlim([0, 8])
ylim([0, 3])



subplot(132)
plot(pfitsall1(:,3), pfitsall1(:,2), 'x');
hold on
plot(pfitsall2(:,3), pfitsall2(:,2), 'x');
plot(pfitsall3(:,3), pfitsall3(:,2), 'x');
plot(pfitsall4(:,3), pfitsall4(:,2), 'x');
xlim([0, .5])
ylim([0, 3])

subplot(133)
plot(pfitsall1(:,1), pfitsall1(:,3), 'x');
hold on
plot(pfitsall2(:,1), pfitsall2(:,3), 'x');
plot(pfitsall3(:,1), pfitsall3(:,3), 'x');
plot(pfitsall4(:,1), pfitsall4(:,3), 'x');
xlim([0, 3])
ylim([0, .5])

%%
% plot the fitted sigmoids
figure;
subplot(221)
% plot(ypredall(manualflat == 1, :)', 'b');
plot(ypredall(classflat1,:)', 'b');

subplot(222)
% plot(ypredall(manualflat == 2, :)', 'b');
plot(ypredall(classflat2,:)', 'b');


subplot(223)
% plot(ypredall(manualflat == 3, :)', 'b');
plot(ypredall(classflat3,:)', 'b');


subplot(224)
% plot(ypredall(manualflat == 4, :)', 'b');
plot(ypredall(classflat4,:)', 'b');


%% 
% plot the raw waveforms
% plot the fitted sigmoids
figure;
subplot(221)
% plot(aggmeansarr(manualflat == 1, :)', 'b');
plot(aggmeansarr(classflat1, :)', 'b');


subplot(222)
% plot(aggmeansarr(manualflat == 2, :)', 'b');
plot(aggmeansarr(classflat2, :)', 'b');

subplot(223)
% plot(aggmeansarr(manualflat == 3, :)', 'b');
plot(aggmeansarr(classflat3, :)', 'b');

subplot(224)
% plot(aggmeansarr(manualflat == 4, :)', 'b');
plot(aggmeansarr(classflat4, :)', 'b');


%%
% figure;
classes = {classflat1, classflat2, classflat3, classflat4};
figure;
hold on
colors = brewermap(4, 'Set1');
colors = colors([2,1,3,4],:);
colors(4,:) = [0,0,0];
for i = 1:4
%     figure;
%     plot(aggmeansarr(classes{i}, :)', 'k', 'LineWidth', 0.25);
    hold on
    meanline = mean(aggmeansarr(classes{i}, :), 1);
    plot(meanline(1:end-1), 'Color', colors(i,:), 'LineWidth', 1);
    
    filename = sprintf('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/hmmblockFigs/cluster%dprofile.pdf', i);
%     saveas(gcf, filename);
    
end


mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)',...
        'xticks', 0:5:15)


%%
figure;
plot(manualflat == 2, '.');
hold on
plot(classflat2 + 0.1, '.');

%%
figure
plot(aggmeansarr((manualflat' ~= 4) & (classflat4 == 1), :)')

%%
aggnorm = aggmeansarr - mean(aggmeansarr, 1);
V = pca(aggnorm);
Xproj = aggnorm * V;
plot(Xproj(classflat1,1), Xproj(classflat1,2), '.')
hold on
plot(Xproj(classflat2,1), Xproj(classflat2,2), '.')
plot(Xproj(classflat3,1), Xproj(classflat3,2), '.')
plot(Xproj(classflat4,1), Xproj(classflat4,2), '.')



function allmeans = getmeans(obs, zstates)
% visualize trials in state
allmeans = [];
for i = 1:4
%     figure;
    obsfilt = obs(zstates == i-1, :);
    allmeans(i,:) = nanmean(obsfilt, 1);
%     imagesc(obsfilt)
end

end

function L = sigmoidloss(params, xvals, y)
ypred = sigmoid(params, xvals);
L = sum((ypred - y).^2);

end

function y = sigmoid(params, x)
mu = params(1);
sigma = params(2);
lapse = params(3);


y = (1-2*lapse) ./ (1 + exp(-(x-mu)/sigma)) + lapse;

end
