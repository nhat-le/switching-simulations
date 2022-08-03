rootdir = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/hmm_analysis_3p';
files = dir(fullfile(rootdir, 'data_3p/070522b/*hmmblockfit_070522b.mat'));

fprintf('Classifying the blockHMM modes...\n')

paths = pathsetup('hmm3p');
expfitdate = '070522';
opts.rootdir = fullfile(paths.expdatapath, expfitdate);
folders = dir(fullfile(opts.rootdir, ...
    sprintf('*hmmblockfit_*%s.mat', expfitdate)));
opts.filter_blocks_by_lengths = 0;
opts.weighted = 1;
opts.python_assist = 0;
opts.effmethod = 'sim';
opts.savefile = 0;
% opts.model_name = 'decoding_common_010522_withsvmMdl_knn_svm_v10_tsne.mat';
opts.model_name = 'decoding_common_010522_withknnMdl_knn_svm_v10_tsne.mat';
opts.svmmodelpath = paths.decodingmodelpath;


[~, aggmeans_native, aggparams_native] = load_params(folders, opts);
load(fullfile(opts.svmmodelpath, opts.model_name), 'Models');


for i = 24
    opts.mdltype = i;
    opts.mdlid = 5;
    [~, ~, statesFlat, ~] = apply_model(aggparams_native, aggmeans_native, opts);
end

% how many states per animal?
paths = pathsetup('matchingsim');

% averaging two seeds
animalModeInfo = struct;
% animal names and K
animalModeInfo.animals = {};
animalModeInfo.K = [];
for i = 1:numel(folders)
    parts = strsplit(folders(i).name, '_');
    animalModeInfo.animals{i} = parts{1};
    K_from_aggparams = size(aggparams_native{i}, 2);
    animalModeInfo.K(i) = K_from_aggparams;
end

assert(sum(animalModeInfo.K) == numel(statesFlat))
% for 070522
statesFlat(1) = 2; %e46: mode 1
statesFlat(14) = 1; %e56: mode 2
statesFlat(20) = 2; %f01: mode 2


% for 070522b
% statesFlat(1) = 2; %e46: mode 1
% statesFlat(18) = 1; %e56: mode 2
% statesFlat(24) = 2; %f01: mode 2


%%

transcurveInfo = struct;

for i = 1:numel(animalModeInfo.K)
    % load file
    K = animalModeInfo.K(i);
    version = '022322_Murat';
    opts.rootdir = fullfile(paths.blockhmmfitpath, version);


    if i == 1
        idstart = 1;
    else
        idstart = sum(animalModeInfo.K(1:i - 1)) + 1;
    end
    idend = sum(animalModeInfo.K(1:i));

    states_animal = statesFlat(idstart:idend);
%     disp(zclassified');
       
    % Load the individual animal blockfit file
    selection = contains({folders.name}, animalModeInfo.animals{i});
    
    assert(sum(selection) == 1) %make sure one and only one animal blockfit file found
    
    load(fullfile(folders(selection).folder, folders(selection).name), 'zstates', 'params', 'lengths', 'transmat', ...
        'block_lens', 'obs');
    parts = strsplit(folders(selection).name, '_');
    animal_name = parts{1};
    transcurveInfo(i).animal = animal_name;

    figure('Name', animal_name);

    for z = 0:max(zstates)
        nexttile
        obs_state = obs(zstates == z, :);
        transcurveInfo(i).transcurve{z+1} = nanmean(obs_state, 1);
        transcurveInfo(i).transclass(z+1) = states_animal(z + 1);
        transcurveInfo(i).params{z+1} = params(:, z+1);
        transcurveInfo(i).obs{z+1} = obs_state;
        plot(transcurveInfo(i).transcurve{z+1})
        title(states_animal(z + 1));


    end


end

%% post-processing
% load("transcurveInfo_070522b.mat")
animals1 = {};
curves1 = [];
fits1 = [];

animals4 = {};
curves4 = [];
fits4 = [];

animals6 = {};
curves6 = [];
fits6 = [];


% just plotting raw data
for i = 1:numel(transcurveInfo)
    for j = 1:numel(transcurveInfo(i).transcurve)
        curve = transcurveInfo(i).transcurve{j};
        params = transcurveInfo(i).params{j};
        curvefit = mathfuncs.sigmoid(0:14, params(1), params(2), params(3));
        state = transcurveInfo(i).transclass(j);
        animal = transcurveInfo(i).animal;

        if state == 1
            animals1{end+1} = animal;
            curves1(end+1, :) = curve;
            fits1(end+1, :) = curvefit;
        elseif state == 4  || state == 3 || state == 2
            animals4{end+1} = animal;
            curves4(end+1, :) = curve;
            fits4(end+1, :) = curvefit;
        elseif state == 5 || state == 6
            animals6{end+1} = animal;
            curves6(end+1, :) = curve;
            fits6(end+1, :) = curvefit;
        else
            error('invalid state')
        end
    end

end

% visualize
figure('Position', [440,461,748,337]);
% subplot(131)
% plot(curves1', 'k')
% 
% subplot(132)
% plot(curves4', 'k')
% 
% subplot(133)
% plot(curves6', 'k')

subplot(131)
plot(fits1', 'k')
hold on
plot(mean(fits1, 1), 'r', 'LineWidth', 2)
ylim([0, 1])
mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)', ...
    'xytitle', 'State 1', 'font_size', 25, 'xticks', 0:5:15)

subplot(132)
plot(fits4', 'k')
hold on
plot(mean(fits4, 1), 'r', 'LineWidth', 2)
ylim([0, 1])

mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)',...
    'xytitle', 'State 2-4', 'font_size', 25, 'xticks', 0:5:15)


subplot(133)
plot(fits6', 'k')
hold on
plot(mean(fits6, 1), 'r', 'LineWidth', 2)
ylim([0, 1])

mymakeaxis('x_label', 'Trials from block start', 'y_label', 'P(Correct)',...
    'xytitle', 'State 5-6', 'font_size', 25, 'xticks', 0:5:15)

%% Plot the HMM identities per animal

% for 070522b
% hmmidentities = [1 2 4 6 0 0; %e46
%     1 1 2 2 6 6; %e56
%     1 2 5 0 0 0; %f02
%     1 1 4 4 4 5; %e53
%     1 1 2 4 6 0; %fh02
%     1 1 2 2 5 6; %fh03
%     1 1 2 2 4 6; %e54
%     1 1 2 4 4 6]; %f01

% for 070522
hmmidentities = [1 2 4 6 0 0; %e46
    1 1 2 2 6 6; %e56
    1 2 5 0 0 0; %f02
    1 4 0 0 0 0; %e53
    1 1 2 4 6 0; %fh02
    1 1 2 2 5 6; %fh03
    1 1 2 2 4 6; %e54
    1 1 2 4 4 6]; %f01

namelst = {'e46', 'e56', 'f02', 'e53', 'fh02', 'fh03', 'e54', 'f01'};


% cmap = brewermap(6, 'Set1');
% cmap = getPyPlot_cMap('spectral', 10);
% cmap = cmap([6,7,8,9,10],:);
% cmap = [0 0 0; cmap];
cmap = [0 0 0;
    86, 180, 233;
    230, 159, 0;
    230, 159, 0;
    230, 159, 0;
    0, 158, 115;
    0, 158, 115] / 255;

figure;
imagesc(hmmidentities);
colormap(cmap);
hold on
vline((0:6) + 0.5, 'k');
hline((0:numel(namelst)) + 0.5, 'k');

axis xy
mymakeaxis('x_label', 'HMM mode', 'y_label', 'Animal', ...
    'xytitle', 'My title', ...
    'font_size', 22, 'xticks', 1:6, 'yticks', 1:numel(namelst), 'yticklabels', namelst)












