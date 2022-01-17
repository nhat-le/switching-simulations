function [out, opts] = load_and_run_tsne(prob, opts_in)
paths = pathsetup('matchingsim');
filepath = fullfile(paths.svmconfigpath, opts_in.version, ...
    sprintf('opts_prob%.1f-*-tsne_final.mat', prob));
files = dir(filepath);
assert(numel(files) == 1);


load(fullfile(files(1).folder, files(1).name), 'opts', 'out');


end


