% for combining the raw data information with hmm state:
% data information given in GLM_clustering/Data folder
% HMM state information given in
% hmm_analysis_3p/data_3p/version/opto_hmm_info.mat file

% load hmm info information
version = '070522';
load(fullfile('data_3p', version, 'opto_hmm_info.mat'));


% load raw data file
data_folder = '/Users/minhnhatle/Documents/ExternalCode/GLM_clustering/Data';
expfolder = '2020-11-21_1_E46';
matfiles = dir(fullfile(data_folder, expfolder, '*Block.mat'));
assert(numel(matfiles) == 1);
load(fullfile(matfiles(1).folder, matfiles(1).name));

% find the corresponding struct in the animalinfo record
parts = strsplit(expfolder, '_');
datestr = parts{1};
sess_order = parts{2};
animal = parts{3};
full_sessname = sprintf('%s_1_%s_Block.mat', datestr, lower(animal));
animalid = find(strcmpi({animalinfo.animal}, animal));
all_sessnames = animalinfo(animalid).sessnames;

for id = 1:size(all_sessnames, 1)
    if strcmpi(all_sessnames(id,:), full_sessname)
        fprintf('sessfound: %d, %s\n', id, all_sessnames(id,:));
        sessid = id;
        break;
    end
end

zstates = animalinfo(animalid).zclassified{sessid};
hmmblocklens = animalinfo(animalid).block_lens{sessid};


% find the block lengths
targets = block.events.contrastLeftValues;
blockstarts = find(diff(targets));
blockstarts = [1 blockstarts + 1 numel(targets) + 1];
blens = diff(blockstarts);

%%
version = '050522';
data_folder = '/Users/minhnhatle/Documents/ExternalCode/GLM_clustering/Data';

folders = dir(data_folder);
for i = 1:numel(folders)
    if ~contains(folders(i).name, '.')
        expfolder = folders(i).name;
        [hmmblocklens, raw_blens, sess_order] = get_block_lengths(expfolder, version);
        parts = strsplit(expfolder, '_');
        rigboxfolder = sprintf('/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/%s/%s', ...
            parts{3}, parts{1});
        recons_blens = reconstruct_block_lengths(rigboxfolder);

        if numel(hmmblocklens) ~= numel(recons_blens)
            fprintf('mismatch: %d - %s\n', i, expfolder)
        end


        check = compare_block_counts(raw_blens, hmmblocklens, sess_order);
        if ~check
            fprintf('%d - %s: failed\n', i, expfolder)
        end
    end
end


%%
expfolder = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/E46/2020-12-08';
blens = reconstruct_block_lengths(expfolder);



function [hmmblocklens, raw_blens, sess_order] = get_block_lengths(expfolder, version)
% Load the hmm block lengths (from the hmm_analysis_3p/data_3p/version/opto_hmm_info.mat)
% and raw block lengths (from the GLM_clustering/Data folder)
% Inputs: expfolder: name of the raw data folder in the GLM_clustering/Data
% Outputs: hmmblocklens, raw_blens: arrays of block lengths 

load(fullfile('data_3p', version, 'opto_hmm_info.mat'));


% load raw data file
data_folder = '/Users/minhnhatle/Documents/ExternalCode/GLM_clustering/Data';
% expfolder = '2020-11-21_1_E46';
matfiles = dir(fullfile(data_folder, expfolder, '*Block.mat'));
assert(numel(matfiles) == 1);

try
    assert(strcmpi(matfiles(1).name(1:end-10), expfolder));
catch
    error('failed')
end

load(fullfile(matfiles(1).folder, matfiles(1).name), 'block');

% find the corresponding struct in the animalinfo record
parts = strsplit(expfolder, '_');
datestr = parts{1};
sess_order = parts{2};
animal = parts{3};
full_sessname = sprintf('%s_1_%s_Block.mat', datestr, lower(animal));

% TODO: parse e53 sessions
if strcmpi(animal, 'e53')
    hmmblocklens = [1 1 1 1];
    raw_blens = [1 1 1 1];
    sess_order = 1;
    return;
end

animalid = find(strcmpi({animalinfo.animal}, animal));
assert(numel(animalid) == 1)
all_sessnames = animalinfo(animalid).sessnames;

for id = 1:size(all_sessnames, 1)
    if strcmpi(all_sessnames(id,:), full_sessname)
%         fprintf('sessfound: %d, %s\n', id, all_sessnames(id,:));
        sessid = id;
        break;
    end
end
hmmblocklens = animalinfo(animalid).block_lens{sessid};


% find the block lengths
targets = block.events.contrastLeftValues;
blockstarts = find(diff(targets));
blockstarts = [1 blockstarts + 1 numel(targets) + 1];
raw_blens = diff(blockstarts);

end



function blens = reconstruct_block_lengths(expfolder)
% Inputs: expfolder: path to the rigbox folder of a single day (can consist
% of multiple sub-sessions 1/2/3 etc
% Returns: blen: concatenate the trials in that day and return the 
% array of block lengths
subdirs = dir(fullfile(expfolder, '*/*Block.mat'));
allchoices = [];
alltargets = [];
allfeedbacks = [];
for i = 1:numel(subdirs)
    load(fullfile(subdirs(i).folder, subdirs(i).name));
    choices = block.events.responseValues;
    targets = block.events.contrastLeftValues;
    feedback = block.events.feedbackValues;
    
    % note: might need to take care of numel(opto)
    N = min([numel(choices) numel(targets) numel(feedback)]);
        
    allchoices = [allchoices choices(1:N)];
    alltargets = [alltargets targets(1:N)];
    allfeedbacks = [allfeedbacks feedback(1:N)];
end

blockstarts = find(diff(alltargets));
blockstarts = [1 blockstarts + 1 numel(alltargets) + 1];
blens = diff(blockstarts);



end






function check = compare_block_counts(blen_raw, blen_hmm, sess_id)
% Compare the block lengths from the raw block.mat file (blen_raw)
% and the block lengths returned by the opto_hmm_info.mat file (blen_hmm)
% sess_id: id of the session (1/2/3)
% returns: 1 if the lengths match

assert(numel(blen_raw) >= 1)
assert(numel(blen_hmm) >= 1)
% assert(numel(blen_hmm) >= numel(blen_raw))

if numel(blen_hmm) < numel(blen_raw)
    assert(numel(blen_hmm) == numel(blen_raw) - 1);
    blen_raw = blen_raw(1:end-1);
end


% Remove singleton trailing blocks
if blen_raw(end) == 1
    blen_raw = blen_raw(1:end-1);
end


if ischar(sess_id)
    sess_id = str2double(sess_id);
end

if sess_id == 1
    % first session of the day, blen_raw should match blen_hmm first
    % entries
    for i = 1:numel(blen_raw) - 1
        if blen_raw(i) ~= blen_hmm(i)
            check = 0;
            fprintf('Length mismatch at position %d\n', i);
            return;
        end
    end


    % last entry might be off by 1...
    if blen_raw(end) > blen_hmm(numel(blen_raw)) + 1
        check = 0;
        fprintf('Length mismatch at position %d\n', numel(blen_raw));
    end

    check = 1;
    return


else
    check = 1; % TODO: implement this
end



end




