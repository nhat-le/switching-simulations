%% Load the simulation results:
qfiledir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyQLearningAgent-withCorr-prob0.20to0.80-072321.mat';
load(qfiledir, 'efflist', 'PLoffsetlist', 'PLslopelist', 'LapseL', 'epslst', 'gammalst');

qefflist = efflist;
qoffset = PLoffsetlist;
qslope = PLslopelist;
qlapse = LapseL;

[epsgrid, gammagrid] = meshgrid(epslst, gammalst);
epsflat = epsgrid(:);
gammaflat = gammagrid(:);

[epsid, gammaid] = meshgrid(1:numel(epslst), 1:numel(gammalst));
epsid_flat = epsid(:);
gammaid_flat = gammaid(:);






%%
ibfiledir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-withCorr-prob0.20to0.80-072121.mat';
load(ibfiledir, 'efflist', 'PLoffsetlist', 'PLslopelist', 'LapseL', 'prewlst', 'pswitchlst');

ibefflist = efflist;
iboffset = PLoffsetlist;
ibslope = PLslopelist;
iblapse = LapseL;

[prewgrid, pswgrid] = meshgrid(prewlst, pswitchlst);
prewflat = prewgrid(:);
pswflat = pswgrid(:);

[prewid, pswid] = meshgrid(1:numel(prewlst), 1:numel(pswitchlst));
prewid_flat = prewid(:);
pswid_flat = pswid(:);


%% Flatten and normalize
[qeffnorm, ibeffnorm] = combine_and_normalize(qefflist, ibefflist);
[qoffsetnorm, iboffsetnorm] = combine_and_normalize(qoffset, iboffset);
[qslopenorm, ibslopenorm] = combine_and_normalize(qslope, ibslope);
[qlapsenorm, iblapsenorm] = combine_and_normalize(qlapse, iblapse);



%% Correlate
D = nan(numel(qeffnorm), numel(ibeffnorm));
for i = 1:numel(qeffnorm)
    for j = 1:numel(ibeffnorm)
        qarr = [qeffnorm(i) qoffsetnorm(i) qslopenorm(i) qlapsenorm(i)];
        ibarr = [ibeffnorm(j) iboffsetnorm(j) ibslopenorm(j) iblapsenorm(j)];
        D(i,j) = sum((qarr - ibarr).^2);
        
    end
end

% For each row, find the idx of the column that is closest to 0
[ibmin_sim, ibnearestidx] = min(D, [], 1);

%idx = 150 seems overwhelming. That corresponds to..eps = 0.139, gamma=1.4
ibnearest = reshape(ibnearestidx, size(ibefflist));
ibmin_sim = reshape(ibmin_sim, size(ibefflist));


[qmin_sim, qnearestidx] = min(D, [], 2);
qmin_sim = reshape(qmin_sim, size(qefflist));


%idx = 150 seems overwhelming. That corresponds to..eps = 0.139, gamma=1.4
qnearest = reshape(qnearestidx, size(qefflist));

%%
% Colormaps for the space of inference-based and q-learning agents
C_ibspace = make_2d_colormap(numel(prewlst), numel(pswitchlst));
C_qspace = make_2d_colormap(numel(epslst), numel(gammalst));

% eps/gamma color
epsid_nearest = epsid_flat(ibnearest);
gammaid_nearest = gammaid_flat(ibnearest);
nearest_color_ib = zeros(size(epsid_nearest,1), size(epsid_nearest, 2), 3);
for i = 1:size(epsid_nearest, 1)
    for j = 1:size(gammaid_nearest, 2)
        nearest_color_ib(i,j,:) = C_qspace(gammaid_nearest(i,j), epsid_nearest(i,j),:);
    end
end


% pr/ps color
prewid_nearest = prewid_flat(qnearest);
pswid_nearest = pswid_flat(qnearest);
nearest_color_q = zeros(size(prewid_nearest,1), size(pswid_nearest, 2), 3);
for i = 1:size(prewid_nearest, 1)
    for j = 1:size(prewid_nearest, 2)
        nearest_color_q(i,j,:) = C_ibspace(pswid_nearest(i,j), prewid_nearest(i,j),:);
    end
end



%%
figure;
subplot(231);
produce_heatmap(qmin_sim, epslst, gammalst, ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', ...
    'clim', [0, 20], 'legendname', 'prew', 'newfig', 0, 'font_size', 12, ...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438])

subplot(232)
produce_heatmap(nearest_color_q, epslst, gammalst, ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', ...
    'clim', [0, 20], 'legendname', 'similarity', 'newfig', 0, 'font_size', 12,...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438])

subplot(233)
produce_heatmap(C_ibspace, prewlst, pswitchlst, ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', ...
    'clim', [0, 1], 'legendname', '\epsilon', 'newfig', 0, 'font_size', 12,...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438], 'ytickvalues', 0:0.1:0.4)


subplot(234)
produce_heatmap(ibmin_sim, prewlst, pswitchlst, ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', ...
    'clim', [0, 20], 'legendname', 'similarity', 'newfig', 0, 'font_size', 12,...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438], 'ytickvalues', 0:0.1:0.4)


subplot(235)
produce_heatmap(nearest_color_ib, prewlst, pswitchlst, ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', ...
    'clim', [0, 1.5], 'legendname', '\gamma', 'newfig', 0, 'font_size', 12,...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438],'ytickvalues', 0:0.1:0.4)

subplot(236)

produce_heatmap(C_qspace, epslst, gammalst, ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$', ...
    'clim', [0, 20], 'legendname', 'similarity', 'newfig', 0, 'font_size', 12,...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438])


%% Individual figures;
produce_heatmap(qmin_sim, epslst, gammalst, 'clim', [0 10], 'legendname', 'Similarity to inference-based', ...
    'x_label', '$\epsilon$', 'y_label', '$\gamma$')

%%
produce_heatmap(ibmin_sim, prewlst, pswitchlst, ...
    'x_label', '$P_{rew}$', 'y_label', '$P_{switch}$', ...
    'clim', [0, 20], 'legendname', 'Similarity to Q-learning', 'font_size', 16,...
    'cbar_pos', [0.3499,0.6791,0.0098,0.2438], 'ytickvalues', 0:0.1:0.4)


%%

function [xnorm, ynorm] = combine_and_normalize(xarr, yarr)
xy_combined = [xarr(:); yarr(:)];
xy_combined = (xy_combined - median(xy_combined)) / std(xy_combined);
xnorm = xy_combined(1:numel(xarr));
ynorm = xy_combined(numel(xarr) + 1 : end);
end


function C = make_2d_colormap(xsize, ysize)
R=[1 0;
   1 0];
G=[1 1
   0 0];
B=[0 0
   0 1];
[xgrid, ygrid] = meshgrid([0, 1], [0, 1]);
[xgridnew, ygridnew] = meshgrid(linspace(0, 1, xsize), ...
    linspace(0, 1, ysize));
Rint = interp2(xgrid, ygrid, R, xgridnew, ygridnew);
Gint = interp2(xgrid, ygrid, G, xgridnew, ygridnew);
Bint = interp2(xgrid, ygrid, B, xgridnew, ygridnew);

C = cat(3, Rint, Gint, Bint);

end


