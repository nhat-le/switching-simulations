% script for analysis of opto animals (f26, 27, 29, 32) 
% modified on 4.2.2022

paths = pathsetup('opto');
% expfitdate = '080122';
expfitdate = '090422';
opts.rootdir = fullfile(paths.opto_expdatapath, expfitdate);
folders = dir(fullfile(opts.rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
mdltypes = 1:2;
mdlids = 1:10;
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.prob = 1; %80-20 world
opts.python_assist = 0;
opts.effmethod = 'sim';
opts.savefile = 1;
% opts.model_name = 'decoding_common_010522_withsvmMdl_knn_svm_v10_tsne.mat';

if opts.prob == 1
    opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10_tsne.mat';
else
    opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10b_tsne.mat';

end
opts.svmmodelpath = paths.decodingmodelpath;


[~, aggmeans_native, aggparams_native] = load_params(folders, opts);


%% Decoding into inf-based/qlearning strategies
fprintf('Classifying the blockHMM modes...\n')

load(fullfile(opts.svmmodelpath, opts.model_name), 'Models');

num_state5 = [];
for i = 24
    opts.mdltype = i;
    opts.mdlid = 5;

    if opts.prob == 1
        [~, ~, statesFlat, ~] = apply_model(aggparams_native, aggmeans_native, opts);
    else
        [~, ~, statesFlat, ~] = apply_model_prob(aggparams_native, aggmeans_native, opts);

    end
    fprintf('i=%d, num states 6 = %d\n', i, sum(statesFlat == 6));
end

%% classify and visualize state distributions
paths = pathsetup('matchingsim');

% averaging two seeds
animalModeInfo = struct;
% animal names and K
animalModeInfo.animals = {};
animalModeInfo.K = [];
for i = 1:numel(folders)
    parts = strsplit(folders(i).name, '_');
    animalModeInfo.animals{i} = parts{1};
    blockfitfile = dir(fullfile(sprintf('%s/%s_hmmblockfit*.mat', ...
        folders(1).folder, parts{1})));
    assert(numel(blockfitfile) == 1);
    load(fullfile(blockfitfile(1).folder, blockfitfile(1).name), 'transmat');
    animalModeInfo.K(i) = size(transmat, 1);

end
animalinfo = struct;

for i = 1:numel(animalModeInfo.K)
    % load file
    K = animalModeInfo.K(i);
    version = '022322_Murat';
    opts.rootdir = fullfile(paths.blockhmmfitpath, version);

       
    % Load the individual animal blockfit file
    selection = contains({folders.name}, animalModeInfo.animals{i});
    
    assert(sum(selection) == 1) %make sure one and only one animal blockfit file found
    
    load(fullfile(folders(selection).folder, folders(selection).name), 'zstates', 'params', 'lengths', 'transmat', ...
        'opto', 'sessnames', 'block_lens', 'obs');
    parts = strsplit(folders(selection).name, '_');
    animal_name = parts{1};
    
    
    assert(strcmp(animal_name, animalModeInfo.animals{i}));
    
    
    zclassified = zstates;
    
    n_zstates = max(zstates) + 1;
    assert(n_zstates == K);
    assert(n_zstates == size(params, 2))
    
    % make sure the states are in the right order
    assert(sum(sum(aggparams_native{selection}(1:3,:) ~= params)) == 0);
    
    counter = K * (find(selection) - 1) + 1;
    statesFlat_extracted = statesFlat(counter : counter + n_zstates - 1);
    for istate = 1:n_zstates
        zclassified(zstates == istate - 1) = statesFlat_extracted(istate);
    end
    
    [identityPermuted, stateSortIdx] = sort(statesFlat_extracted);
    
    
    counter = counter + n_zstates;    
    
    % Break down into individual sessions
    assert(sum(lengths) == numel(zclassified));
    zclassified_splits = mat2cell(zclassified, 1, lengths);
    blocklen_splits = mat2cell(block_lens, 1, lengths);
    obs_splits = mat2cell(obs, lengths);
    
    
    animalinfo(i).animal = animal_name;
    animalinfo(i).zclassified = zclassified_splits;
    animalinfo(i).classes = sort(statesFlat_extracted);
    animalinfo(i).sessnames = sessnames;
    animalinfo(i).opto = cellfun(@(x) find_opto_blocks(x), opto, 'UniformOutput', false);
    animalinfo(i).block_lens = blocklen_splits;
    animalinfo(i).obs_splits = obs_splits;
end



%% Plot state evolution profile for each animal
fprintf('Parsing individual composition...\n');
opts.plotting = 1;
f = waitbar(0);
for id = 1:numel(animalinfo)
    waitbar(id / numel(animalinfo), f);
    zclassified = animalinfo(id).zclassified;
    nstates_regimes = 6; %numel(animalinfo(id).classes);
    composition = nan(numel(zclassified), nstates_regimes);
    for i = 1:numel(zclassified)
        z_single_session = zclassified{i};
        for j = 1:nstates_regimes
            composition(i,j) = sum(z_single_session == j) / numel(z_single_session);   
        end
    end
    
    animalinfo(id).composition = composition;
    nsess = size(composition, 1);
    xlimmax = min(nsess, 37);

    if opts.plotting
        % Plot bar graph to show the composition
        figure('Position', [440,423,736,375], 'Name', animalinfo(id).animal);
        h = bar(composition,'stacked');
        cols = paperaesthetics;
        colors = cols.colors;
%         colors = brewermap(6, 'Set1');
%         orders = [2, 1, 5, 4, 3];
        for i = 1:nstates_regimes
            h(i).FaceColor = colors(i,:);
            h(i).ShowBaseLine = 'off';
        end
        ylim([0 1])
        mymakeaxis('x_label', 'Session #', 'y_label', 'Fraction')

        filename = fullfile(paths.figpath, 'hmmblockFigs/compositions_animals/',...
            sprintf('%s_training_evolution_tsne_021322b_Nmodes.pdf', animalinfo(id).animal));
%         if ~exist(filename, 'file')
%             saveas(gcf, filename);
%         end
    end
end
% close all
close(f)

%% Parse the session info
fprintf('Parsing the log table...\n');
logtbl = readtable(fullfile('optodata/', expfitdate, 'exprecord.xlsx'), ...
    'Sheet', 'Experiments (opto)');

for animalID = 1:numel(animalinfo)
    sessnames = animalinfo(animalID).sessnames;
    name = upper(animalModeInfo.animals{animalID});

    subtbl = logtbl(strcmp(logtbl.Animal_formatted, name), :);
    dates_all = subtbl.Date;
    dates_all.Format = 'dd-MM-yyyy';

    areas_all = {}; %visual/frontal/rsc/motor
    power_all = {}; %high/low
    period_all = {}; %choice/outcome

    for i = 1:size(sessnames, 1)
        % find the corresponding row in logtbl
        filename = sessnames(i,:);
        datename = filename(1:10);

        id = find(dates_all == datename);

        if isempty(id)
            areas_all{i} = '';
            power_all{i} = '';
            period_all{i} = '';
        else
            areas_all{i} = subtbl.Area_formatted{id};
            power_all{i} = subtbl.Power_formatted{id};
            period_all{i} = subtbl.Choice_outcome_formatted{id};
        end



    end

    animalinfo(animalID).areas = areas_all;
    animalinfo(animalID).power = power_all;
    animalinfo(animalID).period = period_all;
end

savefilename = fullfile('optodata/', expfitdate, 'opto_hmm_info.mat');

if exist(savefilename, 'file')
    response = questdlg(sprintf('File exists: %s, overwrite?', savefilename));
    switch response
        case 'Yes'
            save(savefilename, 'opts', 'animalModeInfo', 'folders', 'animalinfo');
            fprintf('File overwritten\n');
        case 'No'
            fprintf('Skipped saving...\n');          
        case 'Cancel'
            error('Operation cancelled')
    end 
    
else
    save(savefilename, 'opts', 'animalModeInfo', 'folders', 'animalinfo')
    fprintf('File saved!\n');
end




function opto = find_opto_blocks(arr)
% Given array, return an array with opto block identities on/off
opto = [];
for i = 1:size(arr, 1)
    block = arr(i,:);
    block = block(~isnan(block));
    opto(i) = block(end);
end

end


