## Notes for running the opto experiment analysis

Note that most opto analysis scripts located in codepath: `MachingSimulations/PaperFigures/code/blockhmm/opto_analysis` 

and most opto processed data located in expdatapath:

`/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data/optodata`

Pipeline:

1. Run `codepath/compile_behavioral_session_opto/`

2. Run `codepath/run_multi_hmm_fits_opto.py`

-- START HERE AFTER PRE-PROCESSING --

3. Run `opto_hmmblockstate_analysis.m` for processing the hmmblockfit files and parsing the information about each session (opto power, region, choice/outcome etc).

The information will be saved in a file in the same folder, `opto_hmm_info.mat`.

4. Run `process_animal_info.m`. This will parse and aggregate the information saved in the opto_hmm_info file and save into a `modecounts.mat` file

5. Run `parse_modecounts.m` for plotting mode fraction with/without opto.
