function [out,opts] = load_data(opts)
% Load features of behavioral simulation
paths = pathsetup('matchingsim');

%% parse options
if ~isfield(opts, 'expfolder')
    opts.expfolder = '092321';
end

if strcmp(opts.expfolder, '092321') || ~isfield(opts, 'rootdir')
    opts.rootdir = paths.simdatapath; %'/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/simdata';
end

if strcmp(opts.expfolder, '092321') || ~isfield(opts, 'filestem')
    opts.filestem{1} = 'EGreedyQLearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat';
    opts.filestem{2} = 'EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-%s.mat';
%     opts.filestem{1} = 'EGreedyqlearningAgent-withCorr-doublesigmoid-prob%.2fto%.2f-100221.mat';
%     opts.filestem{2} = 'EGreedyinf-basedAgent-withCorr-doublesigmoid-prob%.2fto%.2f-100221.mat';
end

if ~isfield(opts, 'prob'); opts.prob = 1; end
if ~isfield(opts, 'outliermode'); opts.outliermode = 1; end

if ~isfield(opts, 'datasavepath')
%     opts.savepath = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';
    opts.datasavepath = sprintf('%s/%s', paths.svmconfigpath, opts.expfolder);
end

if ~isfield(opts, 'figsavepath')
%     opts.savepath = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';
    opts.figsavepath = sprintf('%s/watershedFigs/%s', paths.figpath, opts.expfolder);
end

%%
filedir = sprintf(['%s/%s/' opts.filestem{1}], opts.rootdir, opts.expfolder, 1-opts.prob, opts.prob, opts.expfolder);
load(filedir, 'efflist', 'LapseL', 'PLslopelist', 'PLoffsetlist', 'epslst', 'gammalst');

Qeff_flat = reshape(efflist, [], 1);
Qlapse_flat = reshape(LapseL, [], 1);
Qslope_flat = reshape(PLslopelist, [], 1);
Qoffset_flat = reshape(PLoffsetlist, [], 1);

out.Qdim = size(efflist);

filedir = sprintf(['%s/%s/' opts.filestem{2}], opts.rootdir, opts.expfolder, 1-opts.prob, opts.prob, opts.expfolder);
load(filedir, 'efflist', 'LapseL', 'PLslopelist', 'PLoffsetlist', 'prewlst', 'pswitchlst');

IBeff_flat = reshape(efflist, [], 1);
IBlapse_flat = reshape(LapseL, [], 1);
IBslope_flat = reshape(PLslopelist, [], 1);
IBoffset_flat = reshape(PLoffsetlist, [], 1);

% Filter outliers
if opts.clipmode == 1
    IBoffset_flat(IBoffset_flat < -20) = 3; %-20;
    Qoffset_flat(Qoffset_flat < -20) = 3; %-20;
elseif opts.clipmode == 2
    IBoffset_flat(IBoffset_flat < -20) = -20; %-20;
    Qoffset_flat(Qoffset_flat < -20) = -20; %-20;
end


% Parse the outputs
out.IBdim = size(efflist);

% TODO: change this order to be consistent with the decoding
out.features = [IBeff_flat IBlapse_flat IBslope_flat IBoffset_flat;
    Qeff_flat Qlapse_flat Qslope_flat Qoffset_flat];

out.features_norm = (out.features - nanmean(out.features, 1)) ./ nanstd(out.features, [], 1);

out.prewlst = prewlst;
out.pswitchlst = pswitchlst;
out.gammalst = gammalst;
out.epslst = epslst;

end