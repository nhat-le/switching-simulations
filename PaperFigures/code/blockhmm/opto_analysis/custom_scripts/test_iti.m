filepath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/LocalExpData/f38/2022-06-08/1';
files = dir(fullfile(filepath, '*Block.mat'));

assert(numel(files) == 1);


load(fullfile(files(1).folder, files(1).name), 'block');
assert(contains(block.expDef, 'Choice'));


responseTimes = block.events.responseTimes;
feedbackTimes = block.events.feedbackTimes;
trialStartTimes = block.events.newTrialTimes;

plot(trialStartTimes(2:end) - feedbackTimes(1:end-1));


%%
filepath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/LocalExpData/f39/2022-06-15/1';
files = dir(fullfile(filepath, '*Block.mat'));


assert(numel(files) == 1);

load(fullfile(files(1).folder, files(1).name), 'block');
assert(~contains(block.expDef, 'Choice'));

responseTimes2 = block.events.responseTimes;
feedbackTimes2 = block.events.feedbackTimes;
trialStartTimes2 = block.events.newTrialTimes(1:end-1);



figure;
% plot(trialStartTimes(2:end) - feedbackTimes(1:end-1));
hold on
plot(trialStartTimes2(2:end) - feedbackTimes2(1:end-1));
plot(block.events.contrastLeftValues);





