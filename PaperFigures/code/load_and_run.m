function [out, opts] = load_and_run(prob, varargin)

if numel(varargin) == 0
    usepca = 0;
else
    usepca = varargin{1};
end

folder = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/decodeFigs';

switch prob
    case 0
        filename = 'opts_prob0.0-2021-09-25 20.52.mat';
    case 0.1
        filename = 'opts_prob0.1-2021-09-25 21.44.mat';
    case 0.2
        filename = 'opts_prob0.2-2021-09-25 21.57.mat';
    case 0.3
        filename = 'opts_prob0.3-2021-09-25 22.29.mat';
end

load(fullfile(folder, filename));
opts.save = 0;
opts.savefeatures = 0;
opts.usepca = usepca;


% changed 9.30.21 from run_watershed to run_watershed_pca
if usepca
    [idx, out] = run_watershed_pca(opts);
else
    [idx, out] = run_watershed(opts);
end

% Rotation for idx
if opts.usepca
    switch prob
        case 0
            idx = myrotate(idx, [1,2,3,4]);
        case 0.1
%             idx = myrotate(idx, [6, 1, 3]);
%             idx = myrotate(idx, [4 5]);
        case 0.2
%             idx = myrotate(idx, [2, 1, 4]);
%             idx = myrotate(idx, [5, 3]);
        case 0.3
%             idx = myrotate(idx, [4, 1]);
%             idx = myrotate(idx, [2, 5, 3]);
    end
else
    switch prob
        case 0
            idx = myrotate(idx, [3, 2, 4]);
        case 0.1
            idx = myrotate(idx, [6, 1, 3]);
            idx = myrotate(idx, [4 5]);
        case 0.2
            idx = myrotate(idx, [2, 1, 4]);
            idx = myrotate(idx, [5, 3]);
        case 0.3
            idx = myrotate(idx, [4, 1]);
            idx = myrotate(idx, [2, 5, 3]);
    end
end
%[idxQ, idxIB] = reshapeidx(idx, out);
out.idx = idx;

% Project X on PC space
[~,~,V] = svd(out.features_norm);
out.V = V;

end


