%%
rng(123)
% Make three transition functions
params = [0.5, 2.5, 0.02; %inf-based (fast)
    5, 1, 0.1; %slow switch
    4, 0.1, 0.4]; %random

% Make the functions
N = 15; %number of trials in block
xvals = 0:N-1;
yall = [];
for i = 1:3
    mu = params(i, 1);
    sigma = params(i, 2);
    lapse = params(i, 3);
    yvals = sigmoid(xvals, mu, sigma, lapse);
    yall(:,i) = yvals;
end

% Generate some data
zstates = [1, 1, 1, 2, 1, 1, 3, 3];
feedbacks = [];
for i = 1:numel(zstates)
    z = zstates(i);
    feedbacks(i,:) = rand(N,1) < yall(:,z);   
end

% imagesc(choices);
yflat = reshape(feedbacks', 1, []);

%% Plot raw choices
signedfeedback = feedbacks * 2 - 1;
dirs = repmat([-1; 1], [4, 15]);
choices = signedfeedback .* dirs;
choicesflat = reshape(choices', 1, []);

blockstarts = (0:8) * N;

%plot
figure;
tvals = 1:numel(choicesflat);
plot(tvals(yflat == 1), choicesflat(yflat == 1), 'bo', 'MarkerFaceColor', 'b')
hold on
plot(tvals(yflat == 0), choicesflat(yflat == 0), 'rx', 'LineWidth', 1.5)
vline(blockstarts, 'k--')
colors = brewermap(3, 'Set1');
colors = colors([3, 1, 2],:);
for i = 1:numel(zstates)  
    opts={'EdgeColor', 'none',...
      'FaceColor', colors(zstates(i),:), 'FaceAlpha', 0.2};
    fill_between([blockstarts(i), blockstarts(i+1)], [-1 -1], 1, [], opts{:})
end
ylim([-1.5, 1.5])
mymakeaxis('x_label', 'Trials', 'y_label', 'Choice', 'yticks', [-1, 1])


%% Plot the kernel functions
% figure;
for i = 1:3
    figure;
%     subplot(1,3,i)
    plot(xvals, yall(:,i))
    ylim([0, 1])
    mymakeaxis('x_label', 'Trials in block', 'yticks', [0,1], 'xticks', 0:5:10)
end



%% fh02 transition matrix plot
% figure;
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/fh02_hmmblockfit_102121.mat', 'transmat');

imagesc(transmat)
axis xy
colormap gray
% colorbar

mymakeaxis('x_label', 'z_{t+1}', 'y_label', 'z_{t}', 'xticks', 1:4, 'yticks', 1:4)
caxis([0, 0.6]);
b = colorbar;
b.Ticks = 0:0.1:0.6;
b.Label.String = 'Transition Prob.';
b.Label.FontSize = 16;
b.FontSize = 16;
b.FontName = 'helvetica';
b.Label.FontName = 'helvetica';







function y = sigmoid(x, mu, sigma, lapse)
y = lapse + (1-2*lapse) ./ (1 + exp(-(x-mu) * sigma));

end