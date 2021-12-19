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

Note that the `qlearning_figures_prob.m` might be the more updated script to produce the latest figures that are more convenient to use for all the probabilistic worlds.

Current versions of the figures are based on:
`'/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyQLearningAgent-withCorr-prob0.00to1.00-072321.mat'` and `'/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/simdata/EGreedyInferenceBasedAgent-withCorr-prob0.00to1.00-072121.mat'`

Note that: `PLoffsetlist` is negative, `PLslopelist` is positive in the saved .mat files


### For behavioral regime segmentation with the watershed algorithm:
Run the script `PaperFigures/code/characterization/watershed_segmentation.m`

The following are the standard parameters that were used:
```
opts.nbinhist = 30;
opts.imhmin = 3;
opts.kernelsize = 3
```

Important notes:
* The notebook will save the opts config in the folder `processed_data/svm/configs/{expdate}`
There are option files with the form like `opts_prob0.0-2021-09-25 16.05.mat` which store the parameters of the segmentation, as well as .pdf figures that show the behavioral regime segmentation (this is *not* the final version seen in Fig. 3 in the paper), the MDS space (points in the PC space, fig. 4 in the paper)
  
* Figures are saved in the `PaperFigures/figs/watershedFigs/{expdate}` folder. opts configurations are saved in the `processed_data/svm/configs/{expdate}` folder.

* To visualize the segmentation for all probabilities (final version seen in Fig. 3 of the paper figure), run `PaperFigures/code/characterization/watershed_post_segmentation_plots.m`. This will also produce the segmented spaces that are used in the final figure versions.

* The script also plots the average feature characterization (official figure) + performs the decoding analysis on the short data. It calls the helper function `do_decoding`, works on the saved data from `switching_world_classifier.ipynb` (see below), and then save the decoding models to `processed_data/svm/models`

* Details of watershed pca provided in the `code/decoding/run_watershed_pca.m` script, but main idea is that the projected values Y are obtained by Y = features_norm * V, and features_norm = (features - mean(features)) / std(features)

### For running the decoding analysis which requires simulation on the short data (20-block)
Run `switching_world_classifier.ipynb` notebook to run this simulation
This notebook calls the functions in `decoding.py`. Remember to change the value of `rlow` to the right probability of reward (for e.g. `rlow = 0.1` corresponds to a 90-10 world)

Important notes:
* The notebook will call the opts config stored in the folder `decodeFigs/opts`

* Note that there is a mild filtering taking place in `PaperFigures/code/load_data.m` where offset values below -20 are clipped (disabled this feature on 12.11.2021)

Results will be saved in files with the form `svmresults_from_pickle_092221_prob0.00.mat`, etc.

### For decoding analysis
Run `PaperFigures/code/decoding/decoding_after_watershed.m`

* Need to change the `opts.svmdir` and `opts.svm_version` to the file with the correct date. This is the file where the processed data (Python simulations) are stored (usually in the `processed_data/svm` folder)

* Decoding results and models saved in `processed_data/model` folder

* Note that you will need to change the config file in `processed_data/svm/configs/{expdate}` to the `_final.mat` format so that the file will be recognized. This is because there might be multiple opts files that are saved, and the user should select only one of them to be the final file that will be read for this analysis.

### Schematic of model-free and inf-based agents (Fig. 2 with blue circles and red crosses)

Run `PaperFigures/code/schematic/ibmodels_comparison_figure.m` or `qmodelscomparison_figure.m`

For schematic of the switching dynamics (Fig. 1), run `PaperFigures/code/schematic/qlearning_matching_infbased_schematic.m`
The data for the schematics are located in `processed_data/simdata/schematic`

### BlockHMM
For running the model fitting on simulated data, run `blockHMM.ipynb`

For running the exp fitting procedure, run `blockHMM_expfit.ipynb`
Alternatively, for quickly iterating through all animals, can use the code in the module `src/run_multi_hmm_fits.py`
* If we run this script, the blockhmm fit results will be saved in the folder `processed_data/blockhmmfit/{expdate}`

For running the sigmoidal fits of the blockhmm mode (same fitting procedure as the characterization simulations), we run the notebook `blockhmm_classifier.ipynb`, results will be saved in the `processed_data/blockhmmfit/{expdate}` folder.
* To run this notebook, you need to specify the `version` (the exp date folder that the hmmblockfit data is stored), and the number of states. The script checks that the format of the fit data agrees with this number. Results will be saved as `sigmoid_fit_all__}version}.mat`



### Performing the block-HMM state decoding
Run `code/blockhmm/hmmstate_decode_single_animals.m` for a quick assessment and visualizing the evolution of mixture of strategies of single animal (official figure for fh03 was created here, as well as mean mixture evolution plot)

For visualizing the transition dynamics of each HMM mode, run `code/blockhmm/hmmblockstate_visualize.m` (these are the plots of switching dynamics of the 4 HMM modes in the official figure 8)

`code/blockhmm/hmmblockstate_group.m` is an attempt at grouping the HMM modes in an unbiased way (by calling the decoders that were saved in the previous sections of the paper)

