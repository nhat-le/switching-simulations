% For classifying HMM modes into the behavioral strategy class
% (model-free or inference-based, classes 1 to 5, according to the 
% behavioral model that was fit on simulated data)

%% Load the data
pathsetup;
expfitdate = '113021';
rootdir = fullfile(expdatapath, expfitdate);
folders = dir(fullfile(rootdir, ...
    sprintf('*hmmblockfit_%s.mat', expfitdate)));
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.mdltype = 2;
opts.mdlid = 9;
opts.effmethod = 'rawdata';
opts.model_name = 'decoding_common_121021_withsvmMdl_knn_svm.mat';
opts.svmmodelpath = svmmodelpath;

[params, aggmeans, aggparams] = load_params(folders, opts);
[all_aggparams, aggmeans_all, statesFlat, features_flat] = apply_model(aggparams, aggmeans, opts);

%%
figure(2);
clf;

single_aggparams = {};

for i = 1:5
    single_aggmeans = aggmeans_all(statesFlat == i,:);
    single_aggparams{i} = all_aggparams(:, statesFlat == i);
        
    subplot(2,3,i)
    plot(single_aggmeans', 'b');
    hold on
    plot(mean(single_aggmeans, 1), 'k', 'LineWidth', 3);
    
    
end

%%
figure(1)
hold on
colors = brewermap(5, 'Set1');
handles = [];
for i = 1:5
    arr = single_aggparams{i};
    l = plot(arr(1,:), arr(3,:), 'o', 'Color', colors(i,:),...
        'MarkerFaceColor', colors(i,:));
    if numel(l) > 0
        handles(i) = l;
    end
    
end

legend(handles)


%% Deprecated: visualize effect of blocklens and blockcorrects
blensingle = cell2mat(block_lens_all(statesFlat == 3));
bcorrsingle = cell2mat(block_corr_all(statesFlat == 3));

blencell = block_lens_all(statesFlat == 3);
barrcell = block_corr_all(statesFlat == 3);

blen_nums = cellfun(@numel, blencell);
barr_nums = cellfun(@numel, barrcell);

blocktrans = cumsum(blen_nums);


figure(3)
clf
subplot(211)
plot(blensingle)
hold on
plot(bcorrsingle)
vline(blocktrans)

subplot(212)
plot(double(bcorrsingle) ./ double(blensingle));
l = vline(blocktrans);
set(l, 'LineWidth', 2)

%%
figure;
subplot(211)
plot(blencell{2})
hold on
plot(barrcell{2})

subplot(212)
x1 = double(barrcell{2});
x2 = double(blencell{2});
plot(x1 ./ x2);




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


function [params, aggmeans, aggparams] = load_params(folders, opts)
    effmethod = opts.effmethod;
    filter_blocks_by_lengths = opts.filter_blocks_by_lengths;
    weighted = opts.weighted;


    aggmeans = {};
    aggparams = {};
    block_corr_all = {};
    block_lens_all = {};
    for i = 1:numel(folders)
        load(fullfile(folders(i).folder, folders(i).name));

        if i == numel(folders)
            obs(isnan(obs)) = 1;
        end

        % Mean transition function for all trials in a particular z-state
        allmeans = getmeans(obs, zstates);

        % efficiency second try
        effs = [];
        switch effmethod
            case 'rawdata'
                for zid = 1:4
                    block_corr_filt = double(block_corrs(zstates == zid - 1));
                    block_lens_filt = double(block_lens(zstates == zid - 1));

                    if filter_blocks_by_lengths
                        block_corr_filt = block_corr_filt(block_lens_filt > 15 & block_lens_filt < 25);
                        block_lens_filt = block_lens_filt(block_lens_filt > 15 & block_lens_filt < 25);
                    end

                    if weighted
                        effs(zid) = sum(block_corr_filt) / sum(block_lens_filt);
                    else
                        effs(zid) = mean(block_corr_filt ./ block_lens_filt);
                    end

                    block_corr_all{end+1} = block_corr_filt;
                    block_lens_all{end+1} = block_lens_filt;
                end

            case 'boost'
                for zid = 1:4
                    obsfiltered = obs(zstates == zid-1,:);

                % TODO: determine if this 'boosting' can be improved...
                    effs(zid) = sum(obsfiltered(:) == 1) / numel(obsfiltered) / 15*20;
                end

            case 'sim'
                for zid = 1:4
                    paramset = params(:, zid);
                    delta = 0.1;
                    ntrials = 25;
                    transfunc = mathfuncs.sigmoid(0:delta:ntrials, paramset(1), paramset(2), paramset(3));
                    effs(zid) = sum(transfunc) * delta / ntrials;

                end

        end 

        disp(effs)

        params(end + 1, :) = effs;

        aggmeans{i} = allmeans;
        aggparams{i} = params;   

    end
end

function [all_aggparams, aggmeans_all, statesFlat, features_flat] = apply_model(aggparams, aggmeans, opts)
% Load the model 
% load(fullfile(svmmodelpath, 'decoding_common_101421_withknnMdl.mat'), 'Mdls1');
load(fullfile(opts.svmmodelpath, opts.model_name))

all_aggparams = cell2mat(aggparams);
Mdl = Models{opts.mdltype}{opts.mdlid};
lapseFlat = all_aggparams(3,:)';
effFlat = all_aggparams(4,:)';
offsetFlat = all_aggparams(1,:)';
slopesFlat = -all_aggparams(2,:)';

% Note: no normalization since normalization already handled in svm Mdl
features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];

% features_flat(4, 2) = 0.07;
statesFlat = Mdl.predict(features_flat);
% disp(statesFlat');

% Grouping and plotting the transition function by decoded states
aggmeans_all = cell2mat(aggmeans');
end