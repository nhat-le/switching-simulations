function opts = makeparams_tsne(opts)
% Parameters for tsne
switch opts.prob
    case 1
        opts.perplexity = 30;
        opts.nbinhist = 25;
        opts.imhmin = 3;
        opts.kernelsize = 4;
        opts.seed = 5;
        opts.MaxIter = 1200;
        opts.groups = {[3,5]};
        opts.rotations = {[4,1,6,3], [2,5]};
   
    case 0.9
        opts.perplexity = 25;
        opts.nbinhist = 25;
        opts.imhmin = 3;
        opts.kernelsize = 4;
        opts.seed = 3;
        opts.MaxIter = 1200;
        opts.groups = {[1]};
        opts.rotations = {[3,1,6,5,2,4]};
        
        
    case 0.8
        opts.perplexity = 30;
        opts.nbinhist = 25;
        opts.imhmin = 3;
        opts.kernelsize = 4;
        opts.seed = 5;
        opts.MaxIter = 1120;
        opts.groups = {[5,6]};
        opts.rotations = {[4,1,6,3], [5,2]};
        
        
    case 0.7
        opts.perplexity = 30;
        opts.nbinhist = 25;
        opts.imhmin = 3;
        opts.kernelsize = 4;
        opts.seed = 3;
        opts.MaxIter = 1200;
        opts.groups = {[1]};
        opts.rotations = {[3,1,5], [2,4]};


end