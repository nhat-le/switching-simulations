## Instructions for running the analysis and simulation scripts

### Basic setup
0. Install python, anaconda and jupyter notebook if these have not been installed
1. Clone the directory `git clone https://github.com/nhat-le/switching-simulations`
2. Add the folder `MatchingSimulations/PaperFigures/code/helpers` to the MATLAB path
3. In the script `.../code/helpers/pathsetup.m`, modify the following directories:
   * `out.rigboxpath` should point to the directory where the rigbox behavioral datafiles are stored
   * `out.rootpath` should point to the root directory (the full path to the `MatchingSimulations` folder)
4. Install the following Python packages: `smartload, ssm` 

`python -m pip install [package_name]`

### MATLAB setup


### 0. Experimental schematics
* Fig. 1C (illustration of behavior metrics s, alpha, epsilon, E): run `PaperFigures/code/schematic/setup_figure.m`
* Example behavior of model-free and inference-based agent:

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

### Schematic of model-free and inf-based agents (Fig 3 and 4 with blue circles and red crosses)

Run `PaperFigures/code/schematic/ibmodels_comparison_figure.m` or `qmodelscomparison_figure.m`

For schematic of the switching dynamics (Fig. 1), run `PaperFigures/code/schematic/qlearning_matching_infbased_schematic.m`
The data for the schematics are located in `processed_data/simdata/schematic`

### BlockHMM
First, run script `/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/blockhmm/compile_behavioral_sessions.m` to produce _all_sessions files

Then, modify the `fitrange` file, located in `processed_data/expdata/fitranges_122221.mat`. This file indicates the ranges of files that are good to analyze. For example, if `fitrange` is [6,7,8,9], the blockhmm protocol only looks at the sessions 6 to 9 (note these are Python zero-indexed) in the saved `_all_sessions` file

For running the model fitting on simulated data, run `blockHMM.ipynb`

For running the exp fitting procedure, run `blockHMM_expfit.ipynb`
Alternatively, for quickly iterating through all animals, can use the code in the module `src/run_multi_hmm_fits.py`
* If we run this script, the blockhmm fit results will be saved in the folder `processed_data/blockhmmfit/{expdate}`
* Parameters to change:
   * `version` (line 18): version number of the expdata that is loaded
   * `version_save`  (line 19): version number of the file that is saved (in `blockhmmfit`)
   * `fitrangefile` (line 25): path to the fitrange file which specifies the range of files to fit
    * `savefile` (line 66): whether to save the results

For running the sigmoidal fits of the blockhmm mode (same fitting procedure as the characterization simulations), we run the notebook `blockhmm_classifier.ipynb`, results will be saved in the `processed_data/blockhmmfit/{expdate}` folder.
* To run this notebook, you need to specify the `version` (the exp date folder that the hmmblockfit data is stored), and the number of states. The script checks that the format of the fit data agrees with this number. Results will be saved as `sigmoid_fit_all_{version}.mat`



### Performing the block-HMM state decoding
Run `code/blockhmm/hmmstate_decode_single_animals.m` for a quick assessment and visualizing the evolution of mixture of strategies of single animal (official figure for fh03 was created here, as well as mean mixture evolution plot)

For visualizing the transition dynamics of each HMM mode, run `code/blockhmm/hmmblockstate_visualize_122121.m` (these are the plots of switching dynamics of the 4 HMM modes in the official figure 8)

`code/blockhmm/hmmblockstate_group.m` is an attempt at grouping the HMM modes in an unbiased way (by calling the decoders that were saved in the previous sections of the paper)

### Block HMM selection of K by cross-validation
Run `blockHMM_K_selection.ipynb`. Results are saved as a pkl file, for e.g. `blockhmm_validation_120221.pkl`

