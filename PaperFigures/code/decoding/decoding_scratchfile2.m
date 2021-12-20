% Aim: compare the regime boundaries for nblocks = 1000 vs nblocks = 20
paths = pathsetup('matchingsim');
opts = struct;
opts.rootdir = paths.simdatapath;
opts.expfolder = '121021';
opts.nbinhist = 30;
opts.imhmin = 3;
opts.kernelsize = 3;
opts.prob = 1;
opts.save = 1;
opts.seed = 3;
opts.rotations = {[3, 5], [2, 4]};
opts.plotfeatures = 0;
rng(opts.seed);

% form the feature vectors
[out,opts] = load_data(opts);


%%
paths = pathsetup('matchingsim');
% Load raw data
opts1 = struct;
opts1.version = '121021';
[res1, opts] = load_and_run(0, opts1);


%%
opts2.version = '121721';
[res2, opts] = load_and_run(0, opts2);


%% Re-deriving pcas
features_norm = (res2.features - nanmean(res2.features)) ./ nanstd(res2.features, [], 1);
V = res2.V;
V(:,1) = -V(:,1);
Y = features_norm * V;

% obs = [-1.21, 1.21, 0.01, 0.981]; %need to switch order
obs = [0.981, 0.01, 1.21, -1.21];
obs_norm = (obs - nanmean(res2.features)) ./ nanstd(res2.features, [], 1);
obsY = obs_norm * V;


%%
figure;
hold on

for idx1 = 1:4
    for idx2 = 1:4
        if idx1 <= idx2
            continue
        else
            figure;
            hold on
            for i = 1:5
                plot(res2.features_norm(res1.idx == i,idx1), res2.features_norm(res1.idx == i,idx2), '.');
            end
            plot(obs_norm(idx1), obs_norm(idx2), 'rx')

        end
    end
end










