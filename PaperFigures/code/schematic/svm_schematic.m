% Load the SVM simulation results
paths = pathsetup('matchingsim');
% Load raw data
opts = struct;
opts.version = '092321';

[res, opts_out] = load_and_run_tsne(0, opts);

opts.nClasses = max(res.idx);
opts.reps = 20;
opts.method = 'knn';
opts.nNeighbors = 1;
opts.svmdir = paths.svmdatapath;
opts.svm_version = '010522';
opts.save_model = 0;
opts.prob = 1;

load(sprintf('%s/%s/svmresults_from_pickle_%s_prob%.2f.mat', opts.svmdir,...
    opts.svm_version, opts.svm_version, 1-opts.prob));

[idxQ, idxIB] = reshapeidx(res.idx, res);

%% Plot some example transition functions with example group ids
groupids = [2,1,6,4,2,6];
xvals = 0:15;
transfuncs_all = [];
for i = 1:numel(groupids)
    groupid = groupids(i);
    if groupid <= 4 %qlearning
        [x,y] = find(idxQ == groupid);
        elem = randi(numel(x));
        xsingle = x(elem);
        ysingle = y(elem);
        lapse = Qlapse_arr(xsingle, ysingle);
        offset = Qoffset_arr(xsingle, ysingle);
        slope = Qslope_arr(xsingle, ysingle);
        fprintf('lapse = %.2f, offset = %.2f, slope = %.2f\n', lapse, offset, slope);
        transfuncs_all(i,:) = mathfuncs.sigmoid(xvals, -offset, slope, lapse);
        
    else %ib
        [x,y] = find(idxIB  == groupid);
        elem = randi(numel(x));
        xsingle = x(elem);
        ysingle = y(elem);
        lapse = IBlapse_arr(xsingle, ysingle);
        offset = IBoffset_arr(xsingle, ysingle);
        slope = IBslope_arr(xsingle, ysingle);
        fprintf('lapse = %.2f, offset = %.2f, slope = %.2f\n', lapse, offset, slope);
        transfuncs_all(i,:) = mathfuncs.sigmoid(xvals, -offset, slope, lapse);
        
    end
        
    
    
    
end


%% Visualize
figure;
hold on
nrows = ceil(numel(groupids) / 2);
for i = 1:numel(groupids)
    subplot(2,nrows,i)
    plot(0:15, transfuncs_all(i,:), 'k', 'LineWidth', 1)
    ylim([0, 1])
    mymakeaxis('x_label', 'Trials in block', 'y_label', 'P(Correct)');
    
end

