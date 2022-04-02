% For classifying HMM modes into the behavioral strategy class
% (model-free or inference-based, classes 1 to 5, according to the 
% behavioral model that was fit on simulated data)

%% Load the data
paths = pathsetup('matchingsim');
expfitdate = '121821';
rootdir = fullfile(paths.blockhmmfitpath, expfitdate);
folders = dir(fullfile(rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
mdltypes = 1:2;
mdlids = 1:10;
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.python_assist = 0;
opts.effmethod = 'sim';
opts.model_name = 'decoding_common_121721_withsvmMdl_knn_svm_v5.mat';
opts.svmmodelpath = paths.svmmodelpath;


% For loading python-assisted param fitting
sigmoid_file = dir(sprintf('%s/sigmoid_fit_all_%s.mat', rootdir, expfitdate));
load(fullfile(sigmoid_file(1).folder, sigmoid_file(1).name));

folder_names = {folders.name};
folders_ordered = folders;


% Rearange files and folders to have the same order
for i = 1:size(files, 1)
    file_path = strtrim(files(i,:));
    fileparts = strsplit(file_path, '/');
    filename = fileparts{end};

    subparts = strsplit(filename, '_');
    animal = subparts{1};
    
    % reorder the folders
    idx = contains(folder_names, animal);
    assert(sum(idx) == 1)
    
    folders_ordered(i) = folders(idx);
    
    
    
end


[~, aggmeans_double, aggparams_double] = load_params_python_assisted(files, opts, params_all);

[~, aggmeans_native, aggparams_native] = load_params(folders_ordered, opts);


if opts.python_assist
    aggparams = aggparams_double;
    aggmeans = aggmeans_double;
else
    aggparams = aggparams_native;
    aggmeans = aggmeans_native;
end

%%
load(fullfile(opts.svmmodelpath, opts.model_name), 'Models');
%%
num_state5 = [];
f = waitbar(0);
for mdltype = 7 %1:numel(Models)
    waitbar(mdltype / numel(Models), f);
    for mdlid = 1 %1:numel(Models{1})
        opts.mdltype = mdltype;
        opts.mdlid = mdlid;
        [all_aggparams, aggmeans_all, statesFlat, features_flat] = apply_model(aggparams, aggmeans, opts);

        % Plot
        figure('Name', sprintf('mdltype = %d, mdlid = %d', mdltype, mdlid));
%         clf;

        single_aggparams = {};
        num_state5(end+1) = sum(statesFlat == 5);

        %Plot
        for i = 1:5 %5 performance regimes
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
close(f)

disp(num_state5);


%%
numstate5_arr = reshape(num_state5, [], numel(Models));

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


% function [params, aggmeans, aggparams] = load_params_python_assisted(files, opts, params_all)
%     effmethod = opts.effmethod;
%     filter_blocks_by_lengths = opts.filter_blocks_by_lengths;
%     weighted = opts.weighted;
%     
%     paths = pathsetup('matchingsim');
% 
%     aggmeans = {};
%     aggparams = {};
%     block_corr_all = {};
%     block_lens_all = {};
%     for i = 1:size(files, 1)
%         file_path = strtrim(files(i,:));
%         fileparts = strsplit(file_path, '/');
%         filename = fileparts{end};
%         
%         subparts = strsplit(filename, '_');
%         version = subparts{end}(1:end-4);
%       
%         
%         file_path_full = fullfile(paths.blockhmmfitpath, version, filename);
%         
%         load(file_path_full);
% 
%         % Mean transition function for all trials in a particular z-state
%         allmeans = getmeans(obs, zstates);
%         
%         nstates = size(params, 2);
% 
%         % efficiency second try
%         effs = [];
%         switch effmethod
%             case 'rawdata'
%                 for zid = 1:nstates
%                     block_corr_filt = double(block_corrs(zstates == zid - 1));
%                     block_lens_filt = double(block_lens(zstates == zid - 1));
% 
%                     if filter_blocks_by_lengths
%                         block_corr_filt = block_corr_filt(block_lens_filt > 15 & block_lens_filt < 25);
%                         block_lens_filt = block_lens_filt(block_lens_filt > 15 & block_lens_filt < 25);
%                     end
% 
%                     if weighted
%                         effs(zid) = sum(block_corr_filt) / sum(block_lens_filt);
%                     else
%                         effs(zid) = mean(block_corr_filt ./ block_lens_filt);
%                     end
% 
%                     block_corr_all{end+1} = block_corr_filt;
%                     block_lens_all{end+1} = block_lens_filt;
%                 end
% 
%             case 'boost'
%                 for zid = 1:nstates
%                     obsfiltered = obs(zstates == zid-1,:);
% 
%                 % TODO: determine if this 'boosting' can be improved...
%                     effs(zid) = sum(obsfiltered(:) == 1) / numel(obsfiltered) / 15*20;
%                 end
% 
%             case 'sim'
%                 for zid = 1:nstates
%                     paramset = squeeze(params_all(i, zid, :));
%                     paramset(paramset < -20) = -20;
%                     
%                     delta = 0.1;
%                     ntrials = 25;
%                     transfunc = mathfuncs.sigmoid(0:delta:ntrials, -paramset(1), paramset(2), paramset(3));
%                     effs(zid) = sum(transfunc) * delta / ntrials;
%                 end
% 
%         end 
%         
%         params = squeeze(params_all(i,:,:));
%         params(:,1) = -params(:,1);
%         params(:, end+1) = effs;
% 
%         aggmeans{i} = allmeans;
%         aggparams{i} = params';   
% 
%     end
% end
% 
% 
% 
% function [params, aggmeans, aggparams] = load_params(folders, opts)
%     effmethod = opts.effmethod;
%     filter_blocks_by_lengths = opts.filter_blocks_by_lengths;
%     weighted = opts.weighted;
% 
% 
%     aggmeans = {};
%     aggparams = {};
%     block_corr_all = {};
%     block_lens_all = {};
%     for i = 1:numel(folders)
%         load(fullfile(folders(i).folder, folders(i).name));
% 
%         if i == numel(folders)
%             obs(isnan(obs)) = 1;
%         end
% 
%         % Mean transition function for all trials in a particular z-state
%         allmeans = getmeans(obs, zstates);
% 
%         % efficiency second try
%         effs = [];
%         nstates = size(params, 2);
%         switch effmethod
%             case 'rawdata'
%                 for zid = 1:nstates
%                     block_corr_filt = double(block_corrs(zstates == zid - 1));
%                     block_lens_filt = double(block_lens(zstates == zid - 1));
% 
%                     if filter_blocks_by_lengths
%                         block_corr_filt = block_corr_filt(block_lens_filt > 15 & block_lens_filt < 25);
%                         block_lens_filt = block_lens_filt(block_lens_filt > 15 & block_lens_filt < 25);
%                     end
% 
%                     if weighted
%                         effs(zid) = sum(block_corr_filt) / sum(block_lens_filt);
%                     else
%                         effs(zid) = mean(block_corr_filt ./ block_lens_filt);
%                     end
% 
%                     block_corr_all{end+1} = block_corr_filt;
%                     block_lens_all{end+1} = block_lens_filt;
%                 end
% 
%             case 'boost'
%                 for zid = 1:nstates
%                     obsfiltered = obs(zstates == zid-1,:);
% 
%                 % TODO: determine if this 'boosting' can be improved...
%                     effs(zid) = sum(obsfiltered(:) == 1) / numel(obsfiltered) / 15*20;
%                 end
% 
%             case 'sim'
%                 for zid = 1:nstates
%                     paramset = params(:, zid);
%                     delta = 0.1;
%                     ntrials = 25;
%                     transfunc = mathfuncs.sigmoid(0:delta:ntrials, paramset(1), paramset(2), paramset(3));
%                     effs(zid) = sum(transfunc) * delta / ntrials;
% 
%                 end
% 
%         end 
% 
%         params(end + 1, :) = effs;
% 
%         aggmeans{i} = allmeans;
%         aggparams{i} = params;   
% 
%     end
% end
% 
% 
% function allmeans = getmeans(obs, zstates)
% allmeans = [];
% for i = 1:max(zstates) + 1
%     obsfilt = obs(zstates == i-1, :);
%     allmeans(i,:) = nanmean(obsfilt, 1);
% end
% 
% end
% 
% function [all_aggparams, aggmeans_all, statesFlat, features_flat, MCC] = apply_model(aggparams, aggmeans, opts)
% % Load the model 
% load(fullfile(opts.svmmodelpath, opts.model_name))
% 
% all_aggparams = cell2mat(aggparams);
% Mdl = Models{opts.mdltype}{opts.mdlid};
% MCC = MCCs_all{opts.mdltype}(opts.mdlid);
% offsetFlat = all_aggparams(1,:)';
% slopesFlat = all_aggparams(2,:)';
% lapseFlat = all_aggparams(3,:)';
% effFlat = all_aggparams(4,:)';
% 
% % Note: no normalization since normalization already handled in svm Mdl
% features_flat = [-offsetFlat slopesFlat lapseFlat effFlat];
% 
% features_flat(features_flat < -20) = -20; 
% 
% % features_flat(4, 2) = 0.07;
% statesFlat = Mdl.predict(features_flat);
% 
% % Grouping and plotting the transition function by decoded states
% aggmeans_all = cell2mat(aggmeans');
% end