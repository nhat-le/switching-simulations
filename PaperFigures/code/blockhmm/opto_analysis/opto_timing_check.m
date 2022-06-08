% folder = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/f29/2022-05-07/1';
folder = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/E51/2022-05-09/1';
matfiles = dir(fullfile(folder, '*Block.mat'));
load(fullfile(matfiles(1).folder, matfiles(1).name));

stimOnTimes = block.events.stimulusOnTimes;
feedbackTimes = block.events.feedbackTimes;
targetVals = block.events.contrastLeftValues;
optoVals = block.events.optoblockValues;
N = min([numel(stimOnTimes) numel(feedbackTimes)]);
stimOnTimes = stimOnTimes(1:N);
feedbackTimes = feedbackTimes(1:N);
targetVals = targetVals(1:N);
optoVals = optoVals(2:N);


dt = stimOnTimes(2:end) - feedbackTimes(1:end-1);
dt_fb = feedbackTimes - stimOnTimes;
dt_fb2 = feedbackTimes(2:end) - feedbackTimes(1:end-1);


% block switches
% bswitches = 

%%
figure;
plot(dt)
hold on
plot(targetVals + 1, '.')