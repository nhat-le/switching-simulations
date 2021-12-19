% For classifying HMM modes into the behavioral strategy class
% (model-free or inference-based, classes 1 to 5, according to the 
% behavioral model that was fit on simulated data)

%% Load the data
paths = pathsetup('matchingsim');
expfitdate = '113021';
rootdir = fullfile(paths.blockhmmfitpath, expfitdate);
folders = dir(fullfile(rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
mdltypes = 1:2;
mdlids = 1:10;
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.python_assist = 0;
opts.effmethod = 'sim';
opts.model_name = 'decoding_common_121721_withsvmMdl_knn_svm_v4.mat';
opts.svmmodelpath = paths.svmmodelpath;

if opts.python_assist
    % For loading python-assisted param fitting
    sigmoid_file = dir(fullfile(rootdir, 'sigmoid_fit_all_113021.mat'));
    load(fullfile(sigmoid_file(1).folder, sigmoid_file(1).name));
    [params, aggmeans, aggparams] = load_params_python_assisted(files, opts, params_all);
else
    [params, aggmeans, aggparams] = load_params(folders, opts);
end





%%
all_params = cell2mat(aggparams);
all_params(all_params < -20) = -20;

load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/svm/models/decoding_common_121721_withsvmMdl_knn_svm_v4.mat');
Mdl = Models{1}{2};




%%
% f = waitbar(0);
for mdltype = 1:9 %mdltypes
    for mdlid = 1
        opts.mdltype = mdltype;
        opts.mdlid = mdlid;
%         [all_aggparams, aggmeans_all, statesFlat, features_flat] = apply_model_python_assisted(aggparams, aggmeans, opts);
        [all_aggparams, aggmeans_all, statesFlat, features_flat] = apply_model(aggparams, aggmeans, opts);

        % Plot
        figure('Name', sprintf('mdltype = %d, mdlid = %d', mdltype, mdlid));
        clf;

        single_aggparams = {};
        disp(sum(statesFlat == 5));

        for i = 1:5
            if sum(statesFlat == i) == 0
                continue
            end
            single_aggmeans = aggmeans_all(statesFlat == i,:);
            filtered_aggparams = all_aggparams(:, statesFlat == i);
            filtered_aggparams(filtered_aggparams < -20) = -20;
            single_aggparams{i} = filtered_aggparams;
            subplot(2,3,i)
            plot(single_aggmeans', 'b');
            hold on
            plot(mean(single_aggmeans, 1), 'k', 'LineWidth', 3);
        end
    end
end

%%
% for i =61:80
%     subplot(4,5,i-60)
%     plot(aggmeans_all(i,:))
%     ylim([0, 1])
%     
% end
    
manualclass = [1, 1, 2, 3, 2,...
    1, 1, 3, 2, 1,...
    3, 1, 1, 2, 4,...
    1, 3, 1, 2, 3,...
    2,1,1,3,4,...
    2,3,1,4,2,...
    3,1,3,1,1,...
    2,3,2,2,1,...
    1,3,1,2,3,...
    1,4,2,1,3,...
    2,1,4,2,1,...
    1,3,1,2,4,...
    2,1,1,3,3,...
    1,1,1,1,2,...
    3,1,3,1,4,...
    1];
    
%%
figure;
subplot(221)
plot(aggmeans_all(manualclass==1,:)', 'b')
hold on
plot(mean(aggmeans_all(manualclass==1,:)', 2), 'k', 'LineWidth', 2)

subplot(222)
plot(aggmeans_all(manualclass==2,:)', 'b')
hold on
plot(mean(aggmeans_all(manualclass==2,:)', 2), 'k', 'LineWidth', 2)


subplot(223)
plot(aggmeans_all(manualclass==3,:)', 'b')
hold on
plot(mean(aggmeans_all(manualclass==3,:)', 2), 'k', 'LineWidth', 2)


subplot(224)
plot(aggmeans_all(manualclass==4,:)', 'b')
hold on
plot(mean(aggmeans_all(manualclass==4,:)', 2), 'k', 'LineWidth', 2)


%% Summary of the HMM modes
composition = [];
for i = 1:4
    for j = 1:3
        composition(i,j) = sum(statesFlat == j & manualclass' == i);
    end  
end

statefrac = composition ./ sum(composition, 2);


%%
params_clust3 = features_flat(manualclass == 3,:);


%%
colors = brewermap(6, 'Set1');

figure('Position', [440,423,736,375]);
h = bar(statefrac,'stacked');
orders = [2, 1, 5, 4, 3];
for i = 1:3
    h(i).FaceColor = colors(orders(i),:);
    h(i).ShowBaseLine = 'off';
end
xlim([0.5 4.5])
mymakeaxis('x_label', 'HMM mode', 'y_label', 'Fraction', 'xticks', 1:4)
l = legend(h(1:3), {'Q-learning 1', 'Q-learning 2', 'Q-learning 3'}, 'Position', [0.4,0.42,0.1,0.1], ...
    'FontSize', 14);


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
for i = 1:max(zstates) + 1
%     figure;
    obsfilt = obs(zstates == i-1, :);
    allmeans(i,:) = nanmean(obsfilt, 1);
%     imagesc(obsfilt)
end

end


function [params, aggmeans, aggparams] = load_params_python_assisted(files, opts, params_all)
    effmethod = opts.effmethod;
    filter_blocks_by_lengths = opts.filter_blocks_by_lengths;
    weighted = opts.weighted;
    
    paths = pathsetup('matchingsim');

    aggmeans = {};
    aggparams = {};
    block_corr_all = {};
    block_lens_all = {};
    for i = 1:size(files, 1)
        file_path = strtrim(files(i,:));
        fileparts = strsplit(file_path, '/');
        filename = fileparts{end};
        
        subparts = strsplit(filename, '_');
        version = subparts{end}(1:end-4);
      
        
        file_path_full = fullfile(paths.blockhmmfitpath, version, filename);
        
        load(file_path_full);

        % Mean transition function for all trials in a particular z-state
        allmeans = getmeans(obs, zstates);

        % efficiency second try
        effs = [];
%         effs2 = [];
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
                    paramset = squeeze(params_all(i, zid, :));
                    paramset(paramset < -20) = -20;
                    
%                     paramset1 = paramset(1:3);
%                     paramset2 = paramset(4:6);
                    
                    delta = 0.1;
                    ntrials = 25;
                    transfunc = mathfuncs.sigmoid(0:delta:ntrials, -paramset(1), paramset(2), paramset(3));
                    effs(zid) = sum(transfunc) * delta / ntrials;
%                     effs1(zid) = sum(transfunc1) * delta / ntrials;
%                     transfunc1 = mathfuncs.sigmoid(0:delta:ntrials, -paramset1(1), paramset1(2), paramset1(3));
%                     transfunc2 = mathfuncs.sigmoid(0:delta:ntrials, -paramset2(1), paramset2(2), paramset2(3));
%                     
%                     effs1(zid) = sum(transfunc1) * delta / ntrials;
%                     effs2(zid) = sum(transfunc2) * delta / ntrials;

                end

        end 

%         disp(effs)
        
        params = squeeze(params_all(i,:,:));
%         params1 = params(:,1:3);
%         params2 = params(:,4:6);
%         params = [params1; params2];
        params(:, end+1) = effs;

        aggmeans{i} = allmeans;
        aggparams{i} = params;   

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
        nstates = size(params, 2);
        switch effmethod
            case 'rawdata'
                for zid = 1:nstates
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
                for zid = 1:nstates
                    obsfiltered = obs(zstates == zid-1,:);

                % TODO: determine if this 'boosting' can be improved...
                    effs(zid) = sum(obsfiltered(:) == 1) / numel(obsfiltered) / 15*20;
                end

            case 'sim'
                for zid = 1:nstates
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


function [all_aggparams, aggmeans_all, statesFlat, features_flat, MCC] = apply_model_python_assisted(aggparams, aggmeans, opts)
% Load the model 
% load(fullfile(svmmodelpath, 'decoding_common_101421_withknnMdl.mat'), 'Mdls1');
load(fullfile(opts.svmmodelpath, opts.model_name))

all_aggparams = cell2mat(aggparams');
Mdl = Models{opts.mdltype}{opts.mdlid};
MCC = MCCs_all{opts.mdltype}(opts.mdlid);
offsetFlat = all_aggparams(:,1);
offsetFlat(offsetFlat < -20) = -20; %important!

slopesFlat = all_aggparams(:,2);
lapseFlat = all_aggparams(:,3);
effFlat = all_aggparams(:,4);

% Note: no normalization since normalization already handled in svm Mdl
features_flat = [offsetFlat slopesFlat lapseFlat effFlat];

% features_flat(4, 2) = 0.07;
statesFlat = Mdl.predict(features_flat);
% disp(statesFlat');

% Grouping and plotting the transition function by decoded states
aggmeans_all = cell2mat(aggmeans');
end



function [all_aggparams, aggmeans_all, statesFlat, features_flat, MCC] = apply_model(aggparams, aggmeans, opts)
% Load the model 
% load(fullfile(svmmodelpath, 'decoding_common_101421_withknnMdl.mat'), 'Mdls1');
load(fullfile(opts.svmmodelpath, opts.model_name))

all_aggparams = cell2mat(aggparams);
Mdl = Models{opts.mdltype}{opts.mdlid};
MCC = MCCs_all{opts.mdltype}(opts.mdlid);
offsetFlat = all_aggparams(1,:)';
slopesFlat = all_aggparams(2,:)';
lapseFlat = all_aggparams(3,:)';
effFlat = all_aggparams(4,:)';

% Note: no normalization since normalization already handled in svm Mdl
features_flat = [-offsetFlat slopesFlat lapseFlat effFlat];

% features_flat(4, 2) = 0.07;
statesFlat = Mdl.predict(features_flat);
% disp(statesFlat');

% Grouping and plotting the transition function by decoded states
aggmeans_all = cell2mat(aggmeans');
end