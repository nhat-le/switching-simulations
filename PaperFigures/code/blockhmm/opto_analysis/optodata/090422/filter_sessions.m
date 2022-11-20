animals = {'f32','f33', 'f34', 'f35', 'f36', 'f37', 'f38'};
ranges = {};
tbl = readtable('fullfield_opto_090422.xlsx');
tbl.Date.Format = 'yyyy-MM-dd';
for i = 1:numel(animals)
    fprintf('Processing animal %d of %d: %s...\n', i, numel(animals), animals{i});
    animal = animals{i};
    fname = sprintf('%s_all_sessions_090422.mat', animal);
    load(fname, 'session_names');
    
    idx = filter_animal(animal, tbl);
    ranges{i} = idx - 1;
end


function idx = filter_animal(animal, tbl)
% input: animal: animal name
% tbl: table containing the experimental logs from xls
% output: idx, indices of good sessions to filter

fname = sprintf('%s_all_sessions_090422.mat', animal);
load(fname, 'session_names');
load('fitranges_default.mat', 'animals', 'ranges');

idxanimal = find(strcmpi(animals, animal));
animal_range = ranges{idxanimal};


rigboxfolder = sprintf('/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/%s', animal);

idx = [];

for i = animal_range(1) + 1:numel(session_names)
%     disp(i)
    sess = session_names{i};
    parts = strsplit(sess(1:10), '-');

    % is session 100-0?
    rigboxfilename = dir(fullfile(rigboxfolder, sprintf('*/*/%s', sess)));
    assert(numel(rigboxfilename) == 1)
    try
        load(fullfile(rigboxfilename(1).folder, rigboxfilename(1).name), 'block')
    catch ME
        if (strcmp(ME.identifier,'MATLAB:load:notBinaryFile')) || ...
                (strcmp(ME.identifier,'MATLAB:load:unableToReadMatFile'))
            fprintf('Error loading file: %s\n', rigboxfilename(1).name);
        end  
        continue
    end
    
    % is session in the excel log?
    logrow = find(tbl.Date.Year == str2double(parts{1}) & tbl.Date.Month == str2double(parts{2}) & ...
        tbl.Date.Day == str2double(parts{3}) & strcmpi(tbl.AnimalID, animal));
    if numel(logrow) == 0
        continue
    end


    if contains(block.expDef, 'Prob')
        fprintf('%d: 8020 detected: %s\n', i, sess)
        continue
    end
    
    idx(end+1) = i;
end


% if f32, include also the first sessions (not included in the excel
% spreadsheet)
if strcmpi(animal, 'f32')
    idx = [7:idx(1) - 1, idx];
end



end