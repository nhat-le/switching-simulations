
% animals = {'f01', 'f02', 'f03', 'f04', 'f11', 'f12', 'E35', 'E40',...
%     'fh01', 'fh02', 'f05', 'e53', 'fh03', 'f16', 'f17', 'f20', 'f21', 'f22', 'f23'};
% % animals = {'f25'};
animals = {'f29', 'f26', 'f27', 'f32'};
version = '122221';
paths = pathsetup('matchingsim');
opts.savefile = 1;


f = waitbar(0);
for i = 1:numel(animals)
    waitbar(i/numel(animals), f, sprintf('Processing animal %s', animals{i}));
    opts.root = fullfile(paths.rigboxpath, animals{i});
    opts.savepath = fullfile(paths.expdatapath, version, sprintf('%s_all_sessions_%s.mat', animals{i}, version));
    process_animal(animals{i}, opts); 
end
close(f);


function process_animal(animal, opts)
fprintf('****** Processing animal %s...*******\n', animal);
paths = pathsetup('matchingsim');
root = opts.root;
folders = dir(fullfile(root, '202*'));
choices_cell = {};
targets_cell = {};
feedbacks_cell = {};
maxdelays = [];
probflags = [];
session_names = {};

%TODO: add file names to know experiment type
for id = 1:numel(folders)
%     disp(id)
    [allchoices, alltargets, ~, probflag, allfeedbacks] = rbox.get_choice_sequence(id, root, folders);
    probflags(id) = probflag;
    choices_cell{end+1} = allchoices;
    targets_cell{end+1} = alltargets;
    feedbacks_cell{end+1} = allfeedbacks;
    
    
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
    save(opts.savepath, 'session_names', 'choices_cell', 'targets_cell', 'maxdelays', 'probflags', 'feedbacks_cell');
    fprintf('File saved!\n');
else
    fprintf('Skipping save, file exists...\n');
end

end



