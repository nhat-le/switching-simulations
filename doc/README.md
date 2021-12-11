## Instructions for creating figures in the paper

### 1. For simulations of model-free and inference-based agents (getting the grid metrics with efficiency, slope, lapse, and offset in fig. 2)
* For Python analysis:

Run `switching_world_characterization.ipynb` notebook
Notably, the notebook calls the functions `run_simulations.run_multiple_agents`
Adjust the cell in this notebook which says 
``agent_type = 'qlearning' ``. The type can be `'qlearning'` or `'inf-based'`

#### Notable parameters:
```
- N_iters = 50
- num_states = 2
- obs_dim = 1
- nblocks = 1000
- sigmoid_window = 30
- ntrials_per_block = [10, 40]
```

Results are saved in the files named `EGreedyqlearningAgent-withCorr-doublesigmoid-prob0.00to1.00-121021.mat`

* For MATLAB analysis:

Run `PaperFigures/code/qlearning_figures.m` or `PaperFigures/code/inf_based_figures.m` for the 2d plots of the behavioral metrics

Current versions of the figures are based on:
`'/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat'` and `'/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-withCorr-prob0.00to1.00-072121.mat'`

Note that: `PLoffsetlist` is negative, `PLslopelist` is positive in the saved .mat files


### For behavioral regime segmentation:
Run the script `watershed_segmentation.m`

The following are the standard parameters that were used:
```
opts.nbinhist = 30;
opts.imhmin = 3;
opts.kernelsize = 3
```

Results are saved in the `PaperFigures/decodeFigs` folder. There are option files with the form like `opts_prob0.0-2021-09-25 16.05.mat` which store the parameters of the segmentation, as well as .pdf figures that show the behavioral regime segmentation (Fig. 3 in the paper), the MDS space (points in the PC space, fig. 4 in the paper)

### 2. For running the decoding analysis
Run `switching_world_classifier.ipynb` notebook
This notebook calls the functions in `decoding.py`. Remember to change the value of `rlow` to the right probability of reward (for e.g. `rlow = 0.1` corresponds to a 90-10 world)

Important notes:
* The notebook will call the opts config stored in the folder `decodeFigs/opts`
Results will be saved in files with the form `svmresults_from_pickle_092221_prob0.00.mat`, etc.




