% script for wf animals (e54, e57, f01, f02, f03, f04, f25) made on
% 3.1.2021

paths = pathsetup('matchingsim');
expfitdate = '022822_wf';
opts.rootdir = fullfile(paths.blockhmmfitpath, expfitdate);
folders = dir(fullfile(opts.rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
mdltypes = 1:2;
mdlids = 1:10;
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.python_assist = 0;
opts.effmethod = 'sim';
opts.savefile = 0;
% opts.model_name = 'decoding_common_010522_withsvmMdl_knn_svm_v10_tsne.mat';
opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10_tsne.mat';
opts.svmmodelpath = paths.svmmodelpath;


[~, aggmeans_native, aggparams_native] = load_params(folders, opts);


%% Decoding into inf-based/qlearning strategies
fprintf('Classifying the blockHMM modes...\n')

load(fullfile(opts.svmmodelpath, opts.model_name), 'Models');

num_state5 = [];
for i = 24
    opts.mdltype = i;
    opts.mdlid = 5;
    [~, ~, statesFlat, ~] = apply_model(aggparams_native, aggmeans_native, opts);

    fprintf('i=%d, num states 6 = %d\n', i, sum(statesFlat == 6));
end

%% classify and visualize state distributions
paths = pathsetup('matchingsim');

% averaging two seeds
animalModeInfo = struct;
animalModeInfo.animals = {'e54',...
    'e57',...
    'f01',...
    'f02',...
    'f03',...
    'f04',...
    'f25'};

animalModeInfo.K = [6 6 6 6 6 6 6];
animalinfo = struct;

for i = 1:numel(animalModeInfo.K)
    % load file
    K = animalModeInfo.K(i);
    version = expfitdate;
    opts.rootdir = fullfile(paths.blockhmmfitpath, version);

       
    % Load the individual animal blockfit file
    selection = contains({folders.name}, animalModeInfo.animals{i});
    
    assert(sum(selection) == 1) %make sure one and only one animal blockfit file found
    
    load(fullfile(folders(selection).folder, folders(selection).name), 'zstates', 'params', 'lengths', 'transmat', ...
        'opto', 'sessnames');
    parts = strsplit(folders(selection).name, '_');
    animal_name = parts{1};
    
    
    assert(strcmp(animal_name, animalModeInfo.animals{i}));
    
    
    zclassified = zstates;
    
    n_zstates = max(zstates) + 1;
    assert(n_zstates == K);
    assert(n_zstates == size(params, 2))
    
    % make sure the states are in the right order
    % replace the zstates (raw HMM id) with their decoded regime
    assert(sum(sum(aggparams_native{selection}(1:3,:) ~= params)) == 0);
    
    counter = K * (find(selection) - 1) + 1;
    statesFlat_extracted = statesFlat(counter : counter + n_zstates - 1);
    for istate = 1:n_zstates
        zclassified(zstates == istate - 1) = statesFlat_extracted(istate);
    end
    
    % Permute the states, sorting by the decoded regime
    [identityPermuted, stateSortIdx] = sort(statesFlat_extracted);
    paramsPermuted = params(:, stateSortIdx);
    
    
    counter = counter + n_zstates;    
    
    % Break down into individual sessions
    assert(sum(lengths) == numel(zclassified));
    zclassified_splits = mat2cell(zclassified, 1, lengths);
    zstates_splits = mat2cell(zstates, 1, lengths);
    
    
    animalinfo(i).animal = animal_name;
    animalinfo(i).zraw = zstates_splits;
    animalinfo(i).params = params;
    animalinfo(i).zclassified = zclassified_splits;
    animalinfo(i).classes = identityPermuted;
    animalinfo(i).classes_raw = statesFlat_extracted;
    animalinfo(i).sessnames = sessnames;
    animalinfo(i).opto = cellfun(@(x) x(:,2), opto, 'UniformOutput', false);
    
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
        if ~exist(filename, 'file')
            saveas(gcf, filename);
        end
    end
end
% close all
close(f)

%% Parse the session info
fprintf('Parsing the log table...\n');
logtbl = readtable('experimental_records_raw.xlsx', 'Sheet', 'Experiments (opto)');

for animalID = 1:4
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

savefilename = 'opto_hmm_info_wf.mat';
if ~exist(savefilename, 'file')
    save(savefilename, 'opts', 'animalModeInfo', 'folders', 'animalinfo')
    fprintf('File saved!\n');
end

