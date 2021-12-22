## Notes on data format 

### Data format of simulation params (1000 blocks)
`res.features` is N x 4 matrix (both model-free and inference-based modes are included)
**CAUTION** order is not the same as the decoding rows!
* Column 1: efficiency (positive)
* Column 2: lapse (positive, from 0 to 0.5)
* Column 3: slope (positive, capped at 3)
* Column 4: offset (negative, no cap)

### Data format of blockHMM params
`params` is a 3x4 matrix (3 x n-HMM blocks (K)):

* First row: offset (positive). **CAUTION** to change to negative if used with decoders
* Second row: slope (positive)
* Third row: lapse (between 0 and 0.5)

`code/blockhmm/hmmblockstate_group.m` appends one row to this with efficiency (positive from 0 to 1)

### Data format of raw features (result of `load_and_run`)
* res.features is a N x 4 matrix
* First column: efficiency (from 0.5 to 0.9)
* Second column: lapse (from 0 to 0.5)
* Third column: slope (from 0 to 3, positive)
* Fourth column: offset (negative, capped at -20)

### Data format of decoding params (result of `decoding_after_watershed`)
`features` is a N x 4 matrix

* First column: offset (negative, capped at -20, SVM performance improves with this cap)
* Second column: slope (positive, capped at 3)
* Third column: lapse (positive, between 0 and 0.5)
* Fourth column: efficiency (positive, between 0.5 and 1)

### Data format of blockHMM aggregate params
This is used in the script `hmmblockstate_group`
- aggparams_double contains cell arrays of 4 x `nstates`

* First row: offset (positive)
* Second row: slope (positive)
* Third row: lapse (0 to 0.5)
* Fourth row: efficiency