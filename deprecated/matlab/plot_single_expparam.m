load expdata/expfit_params.mat

Nsess = 20; %plot only this many sessions
bluecol = [33,113,181]/255;


% plot for f16
id = 5;
figure;
plot(expeff_all(id,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
% Fit the eff mean
xvals = find(~isnan(expeff_all(id,:)));
xvals = xvals(1:end);
yvals = expeff_all(id, xvals);
p = [1, 10, 0.2];
pfit = fminsearch(@(p) losssigmoid(p, xvals, yvals), p);
ypred = 0.5 + (1 - 0.5 - pfit(3)) ./ (1 + exp(-pfit(1) * (xvals - pfit(2))));
hold on
plot(xvals, ypred, 'Color', bluecol, 'LineWidth', 2);

mymakeaxis('x_label', 'Session #', 'y_label', 'Performance', 'font_size', 16)


%%
figure;
plot(explapses_all(id,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')

% Fit the offset mean
xvals = find(~isnan(explapses_all(id,:)));
xvals = xvals(2:end);
yvals = explapses_all(id, xvals);
p = [1, 1, 2];
lossfun(p, xvals, yvals)
pfit = fminsearch(@(p) lossfun(p, xvals, yvals), p);
ypred = pfit(1) * exp(-pfit(2)*xvals) + pfit(3);
hold on
plot(xvals, ypred, 'Color', bluecol, 'LineWidth', 2);

mymakeaxis('x_label', 'Session #', 'y_label', 'Lapse', 'font_size', 16)

%%
figure;
plot(expoffsets_all(id,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
% Fit the offset mean
xvals = find(~isnan(expoffsets_all(id,:)));
xvals = xvals(2:end);
yvals = expoffsets_all(id, xvals);
p = [1, 1, 2];
lossfun(p, xvals, yvals)
pfit = fminsearch(@(p) lossfun(p, xvals, yvals), p);
ypred = pfit(1) * exp(-pfit(2)*xvals) + pfit(3);
hold on
plot(xvals, ypred, 'Color', bluecol, 'LineWidth', 2);

mymakeaxis('x_label', 'Session #', 'y_label', 'Offset', 'font_size', 16)

%%
figure;
plot(expslopes_all(id,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
% Linear regression 
xvals = find(~isnan(expslopes_all(id,:)));
xvals = xvals(1:end);
yvals = expslopes_all(id, xvals);
X = [ones(length(xvals),1) xvals'];
b = X\yvals';
ypred = b(1) + b(2) * xvals;
hold on
plot(xvals, ypred, 'Color', bluecol, 'LineWidth', 2);
mymakeaxis('x_label', 'Session #', 'y_label', 'Slope', 'font_size', 16)

function loss = lossfun(p, x, y)
A = p(1);
k = p(2);
c = p(3);
ypred = A * exp(-k*x) + c;
loss = sum((ypred - y).^2);
end


function loss = losssigmoid(p, x, y)
lapseL = 0.5;
slope = p(1);
offset = p(2);
lapseR = p(3);
ypred = lapseL + (1 - lapseL - lapseR) ./ (1 + exp(-slope * (x - offset)));
loss = sum((ypred - y).^2);

end