% behavFolder = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/E45/2020-12-03/1/';
behavFolder = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/F11/2021-09-14/1/';
files = dir(fullfile(behavFolder, '*.mat'));

load(fullfile(files(1).folder, files(1).name));

figure('Position', [90,335,1325,463]);
ax=axes;
plotPerfProcessedHelper(ax,block)

set(gca, 'FontSize', 16);
xlim([0, 259])

mymakeaxis('x_label', 'Trial #', 'y_label', 'Choice', 'xticks', 0:50:250,...
    'yticks', [-1 1], 'yticklabels', {'Left', 'Right'}, 'font_size', 20)


function plotPerfProcessedHelper(ax, block)
%PLOTPERF Summary of this function goes here
%   Detailed explanation goes here

contrast = [];
contrast(1,:) = [block.events.contrastLeftValues];
contrast(2,:) = [block.events.contrastRightValues];

response = [block.events.responseValues];
repeatNum = [block.events.repeatNumValues];
trialside = [block.events.trialSideValues];


% Make contrast, response and repeat num the same length
N = numel(response);
contrast = contrast(:,1:N);
repeatNum = repeatNum(:,1:N);
incl = ~any(isnan([contrast;response;repeatNum]));
contrast = contrast(:,incl);
response = response(incl);
trialside = trialside(incl);
% repeatNum = repeatNum(incl);

% trialstart = 1:numel(response);
dcontrast = diff(contrast, [], 1);

%filter out
filtrange = 76:335;
dcontrast = dcontrast(:, filtrange);
response = response(filtrange);
trialside = trialside(filtrange);
trialstart = 1:numel(response);


% Find choice types

% First plot the time-outs
cols = paperaesthetics;

timeouts = find(response == 0);
plot(ax, trialstart(timeouts), dcontrast(timeouts), 'ko', 'MarkerFaceColor', 'k');
hold(ax, 'on')

% Plot the leftward/rightward trials
leftCorr = find(response == -1 & trialside == response);
leftIncorr = find(response == -1 & trialside ~= response);
rightCorr = find(response == 1 & trialside == response);
rightIncorr = find(response == 1 & trialside ~= response);

markersize = 10;

plot(ax, trialstart(leftCorr), ones(1,numel(leftCorr)) * -1, 'o',...
    'MarkerEdgeColor', 'w', 'MarkerFaceColor', cols.bluecol, 'MarkerSize', markersize)
plot(ax, trialstart(rightCorr), ones(1,numel(rightCorr)), 'o',...
    'MarkerEdgeColor', 'w', 'MarkerFaceColor', cols.bluecol, 'MarkerSize', markersize)
plot(ax, trialstart(leftIncorr), ones(1,numel(leftIncorr)) * -1, 'x', 'Color', cols.redcol,...
     'MarkerSize', markersize)
plot(ax, trialstart(rightIncorr), ones(1,numel(rightIncorr)), 'x', ...
    'Color', cols.redcol, 'MarkerSize', markersize)


% Plot the moving average
N = 15; % define window
convFilter = ones(1,N) / N;
movAverage = conv(response, convFilter, 'same');
plot(ax, movAverage, 'k', 'LineWidth', 0.5)


% Fill areas to indicate block transitions
tTransitions = find(diff(dcontrast));
tTransitions = [0 tTransitions numel(dcontrast)];
for i = 1:numel(tTransitions)-1
    plot(ax, [tTransitions(i) tTransitions(i)], [-2 1.5], 'k--');
    curr = tTransitions(i);
    next = tTransitions(i+1);
    
    if dcontrast(tTransitions(i) + 1) == 1
        fill(ax, [curr next next curr], [2 2 2.2 2.2], 'k', 'EdgeAlpha', 0);       
%         text(ax, (curr + next)/2 - 3, 2.1, 'R', 'Color', 'w', 'FontSize', 16);
    else
        fill(ax, [curr next next curr], [2 2 2.2 2.2], 'w', 'EdgeAlpha', 0);
    end
    plot(ax, [curr curr], [2 2.2], 'k');
    plot(ax, [next next], [2 2.2], 'k');
end

hline([2 2.2], 'k')
% vline(tTransitions, 'k')
% plot([0 numel(response)], [2 2], 'k')


set(ax, 'YTick', [-1 1])
set(ax, 'YTickLabel', {'Left', 'Right'})
ylim(ax, [-3 3])


end

