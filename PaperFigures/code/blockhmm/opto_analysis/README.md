## Notes for running the opto experiment analysis
1. Run the script `opto_hmmblockstate_analysis.m` for processing the hmmblockfit files and parsing the information about each session (opto power, region, choice/outcome etc).

The information will be saved in a file in the same folder, `opto_hmm_info_{version}.mat`.

2. Run `process_animal_info.m`. This will parse and aggregate the information saved in the opto_hmm_info file and save into a `modecounts.mat` file


