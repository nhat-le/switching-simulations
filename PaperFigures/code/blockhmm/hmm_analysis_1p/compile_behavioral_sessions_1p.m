animals = {'e54', 'f01', 'f02', 'f03', 'f04', 'e57', 'f25'};
opts.version = '050522';
paths = pathsetup('hmm1p');
opts.savefile = 1;

f = waitbar(0);

% copy the fitrange file if not exists
savefolder = fullfile(paths.expdatapath, opts.version);
fitrangepath = sprintf('%s/fitranges_%s.mat', savefolder, opts.version);

if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

if ~exist(fitrangepath, 'file')
    % copy default fitrange file
    fprintf('Failed to find fitrange.mat file, copying default file...\n')
    copyfile(paths.default_fitrange, fitrangepath); 
end

%%
for i = 7 %1 :numel(animals)
    waitbar(i/numel(animals), f, sprintf('Processing animal %s', animals{i}));
    opts.root = fullfile(paths.rigboxpath, animals{i});
    opts.savepath = fullfile(paths.expdatapath, opts.version, sprintf('%s_all_sessions_%s.mat', animals{i}, opts.version));
    process_animal(animals{i}, opts); 
end
close(f);


function process_animal(animal, opts)
fprintf('****** Processing animal %s...*******\n', animal);

% check if there is a fitrange file in the folder
savefolder = fileparts(opts.savepath);
fitrangepath = sprintf('%s/fitranges_%s.mat', savefolder, opts.version);
if ~exist(fitrangepath, 'file')
    error('Failed to find fitrange.mat file')
end


root = opts.root;
folders = dir(fullfile(root, '202*'));
choices_cell = {};
targets_cell = {};
feedbacks_cell = {};
sesscounts_cell = {};
opto_cell = {};
maxdelays = [];
probflags = [];
session_names = {};

%TODO: add file names to know experiment type
for id = 1:numel(folders)
%     disp(id)
    [allchoices, alltargets, ~, probflag, allfeedbacks, allopto, allsesscounts] = rbox.get_choice_sequence_opto(id, root, folders);
    probflags(id) = probflag;
    choices_cell{end+1} = allchoices;
    targets_cell{end+1} = alltargets;
    feedbacks_cell{end+1} = allfeedbacks;
    sesscounts_cell{end+1} = allsesscounts;
    opto_cell{end+1} = allopto;
    
    
    files = dir(fullfile(root, ...
        folders(id).name, '*/*Block.mat'));
    session_names{end+1} = files(1).name;
    
    try
        load(fullfile(files(1).folder, files(1).name));
    
        if isfield(block.paramsValues, 'rewardDelay')
            maxdelays(end+1) = max([block.paramsValues.rewardDelay]);
        else
            maxdelays(end+1) = 0;
        end
    catch ME
        if (strcmp(ME.identifier,'MATLAB:load:notBinaryFile')) || ...
                (strcmp(ME.identifier,'MATLAB:load:unableToReadMatFile'))
            fprintf('Error loading file, skipping,...\n');
        end     
        maxdelays(end+1) = nan;
    continue
    end
    
       
end

%%
if opts.savefile && ~exist(opts.savepath, 'file')
    % Create folder if not exists
    if ~exist(fileparts(opts.savepath), 'dir')
        mkdir(fileparts(opts.savepath))
    end
    
    save(opts.savepath, 'session_names', 'choices_cell', 'targets_cell', 'maxdelays', 'probflags', 'feedbacks_cell',...
        'opto_cell', 'sesscounts_cell');
    fprintf('File saved!\n');
    
    
    % update the fitranges file
    load(fitrangepath, 'ranges', 'animals');
    nsess = numel(choices_cell);
       
    idx = find(contains(animals, lower(animal)));
    assert(numel(idx) == 1);
    if nsess - 1 > ranges{idx}(end)
        curr_range = ranges{idx};
        nextras = nsess - curr_range(end) - 1;
        ranges{idx}(end + 1 : end + nextras) = curr_range(end) + 1 : nsess - 1;
        save(fitrangepath, 'ranges', 'animals');
    end

    
elseif opts.savefile
    fprintf('Skipping save, file exists...\n');
end

end



