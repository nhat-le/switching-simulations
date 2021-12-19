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


%% Check
for i = 1:numel(aggmeans_double)
    diff = (aggmeans_double{i} - aggmeans_native{i}).^2;
    assert(sum(diff(:)) == 0);
end

%% Which is the better fit?
for id = 1:17
    aggmeans = aggmeans_double{id}; %they should be the same as native anyway
    aggparams_double_single = aggparams_double{id};
    aggparams_native_single = aggparams_native{id};

    xvals = 1:size(aggmeans, 2);
    nstates = size(aggmeans, 1);

    figure;
    for clustid = 1:nstates
        double_params = aggparams_double_single(:,clustid);
        native_params = aggparams_native_single(:,clustid);

        y_double = mathfuncs.sigmoid(xvals, double_params(1), double_params(2), double_params(3));
        y_native = mathfuncs.sigmoid(xvals, native_params(1), native_params(2), native_params(3));

        subplot(2,3,clustid)
        plot(aggmeans(clustid,:))
        hold on
        l1 = plot(y_double);
        l2 = plot(y_native);
%         legend([l1, l2], {'Double', 'Native'});
    end

end

%% Comparison of double vs native
for i = 1:numel(aggparams_double)
    figure;
    subplot(121)
    imagesc(aggparams_double{i})
    caxis([0 3])
    
    subplot(122)
    imagesc(aggparams_native{i})
    caxis([0 3])
    
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
        
        nstates = size(params, 2);

        % efficiency second try
        effs = [];
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
                    paramset = squeeze(params_all(i, zid, :));
                    paramset(paramset < -20) = -20;
                    
                    delta = 0.1;
                    ntrials = 25;
                    transfunc = mathfuncs.sigmoid(0:delta:ntrials, -paramset(1), paramset(2), paramset(3));
                    effs(zid) = sum(transfunc) * delta / ntrials;
                end

        end 
        
        params = squeeze(params_all(i,:,:));
        params(:,1) = -params(:,1);
        params(:, end+1) = effs;

        aggmeans{i} = allmeans;
        aggparams{i} = params';   

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

        params(end + 1, :) = effs;

        aggmeans{i} = allmeans;
        aggparams{i} = params;   

    end
end


function allmeans = getmeans(obs, zstates)
allmeans = [];
for i = 1:max(zstates) + 1
    obsfilt = obs(zstates == i-1, :);
    allmeans(i,:) = nanmean(obsfilt, 1);
end

end
