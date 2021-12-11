% For classifying HMM modes into the behavioral strategy class
% (model-free or inference-based, classes 1 to 5, according to the 
% behavioral model that was fit on simulated data)

%% Load the data
rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/';
folders = dir(fullfile(rootdir, '*hmmblockfit_113021.mat'));

effmethod = 'sim';

aggmeans = {};
aggparams = {};
block_corr_all = {};
block_lens_all = {};
for i = 1:numel(folders)
    load(fullfile(rootdir, folders(i).name));
    
    if i == numel(folders)
        obs(isnan(obs)) = 1;
    end
    
    % Mean transition function for all trials in a particular z-state
    allmeans = getmeans(obs, zstates);
    
    % efficiency second try
    effs = [];
%     block_corr_single = {};
%     block_lens_single = {};

    switch effmethod
        case 'rawdata'
            for zid = 1:4
                block_corr_filt = block_corrs(zstates == zid - 1);
                block_lens_filt = block_lens(zstates == zid - 1);
                x1 = double(block_corr_filt);
                x2 = double(block_lens_filt);

                effs(zid) = mean(x1(x2 > 15 & x2 < 25) ./ x2(x2 > 15 & x2 < 25));
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
    
%     block_corr_all{i} = block_corr_single;
%     block_lens_all{i} = block_lens_single;
    % Calculate efficiency
%     effs = [];
%     


%% Decoding of aggparams
% Load the model 
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/decoding_common_101421_withknnMdl.mat', 'Mdls1');

all_aggparams = cell2mat(aggparams);
Mdl = Mdls1{1};
lapseFlat = all_aggparams(3,:)';
effFlat = all_aggparams(4,:)';
offsetFlat = all_aggparams(1,:)';
slopesFlat = -all_aggparams(2,:)';

% Note: no normalization since normalization already handled in svm Mdl
features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];

% features_flat(4, 2) = 0.07;
statesFlat = Mdl.predict(features_flat);
% disp(statesFlat');

%%
figure(4)
plot(all_aggparams(3, statesFlat == 3), all_aggparams(4, statesFlat == 3), 'o')


%% Grouping and plotting the transition function by decoded states
aggmeans_all = cell2mat(aggmeans');

figure(2);
clf;

single_aggparams = {};

for i = 1:5
    single_aggmeans = aggmeans_all(statesFlat == i,:);
    single_aggparams{i} = all_aggparams(:, statesFlat == i);
    
%     for blockid = 1:all_aggparams
    
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


%%
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