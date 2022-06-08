% script for simple exploration of 80-20 behavior sessions
% filepath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/LocalExpData/f29/2022-05-26/1';
filepath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/LocalExpData/f35/2022-05-29/1';

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

% blocktrans = blocktrans(1:end-1);

% Get the transition function
perf_all = [];
for i = 0:window
    perf_all(:,i+1) = responses(blocktrans + i);
end

if targets(1) == 1
    perf_all(1:2:end,:) = perf_all(1:2:end,:) * -1;
else
    perf_all(2:2:end,:) = perf_all(2:2:end,:) * -1;
end

figure;
imagesc(perf_all)
colormap gray
title(sprintf('%s behavior %s', animal, expdate))
set(gca, 'FontSize', 20)
xlabel('Trials in block')
ylabel('Block #')

%% Plot average transition
figure;
plot(mean(perf_all > 0), 'k', 'LineWidth', 2)
title(sprintf('%s behavior %s', animal, expdate))
set(gca, 'FontSize', 20)
hline(0.8, 'k--')
xlabel('Trials in block')
ylabel('P(correct)')










