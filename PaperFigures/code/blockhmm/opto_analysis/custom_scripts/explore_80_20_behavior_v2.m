% script for simple exploration of 80-20 behavior sessions
% same as v1 but now use a 'padding' approach
filepath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/LocalExpData/f27/2022-05-27/1';
% filepath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/LocalExpData/f35/2022-05-29/1';

fileparts = strsplit(filepath, '/');
animal = fileparts{end-2};
expdate = fileparts{end-1};

files = dir(fullfile(filepath, '*Block.mat'));
assert(numel(files) == 1)
load(fullfile(files(1).folder, files(1).name))

targets = block.events.contrastLeftValues;
responses = block.events.responseValues;
feedback = block.events.feedbackValues;

blocktrans = find(diff(targets) ~= 0);
blocktrans = [1 blocktrans + 1];

window = 20;
if numel(responses) - blocktrans(end) + 1 < window
    blocktrans = blocktrans(1:end-1);
end

blocktranswithEnd = [blocktrans numel(responses) + 1];


% Find the longest block
blocklens = diff(blocktranswithEnd);
maxlen = max(blocklens);


% Get the transition function
perf_all = nan(numel(blocktrans), maxlen);
for i = 1:numel(blocktranswithEnd)-1
    perf_all(i, 1:blocklens(i)) = responses(blocktranswithEnd(i) : blocktranswithEnd(i+1) - 1);

end

if targets(1) == 1
    perf_all(1:2:end,:) = perf_all(1:2:end,:) * -1;
else
    perf_all(2:2:end,:) = perf_all(2:2:end,:) * -1;
end

perf_viz = perf_all;
perf_viz(isnan(perf_viz)) = 0;

figure;
imagesc(perf_viz)
colormap gray
title(sprintf('%s behavior %s', animal, expdate))
set(gca, 'FontSize', 20)
xlabel('Trials in block')
ylabel('Block #')

%% Plot average transition
figure;
window = 25;
plot(nanmean((perf_all(:,1:window) + 1) / 2), 'k', 'LineWidth', 2)
title(sprintf('%s behavior %s', animal, expdate))
set(gca, 'FontSize', 20)
hline(0.8, 'k--')
xlabel('Trials in block')
ylabel('P(correct)')










