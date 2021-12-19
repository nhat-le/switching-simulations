import autograd.numpy as np
import glob
import autograd.numpy.random as npr
import os
import scipy.io
import numpy as np2
npr.seed(0)
from src.utils import pathsetup

import ssm
import smartload.smartload as smart
from ssm.exputils import load_multiple_sessions, make_savedict
npr.seed(0)

def run_and_save(animal, seed, N_iters=3000, num_states=6):
    print(f'Starting run and save for {animal}, seed {seed}')
    # Load data
    version = '113021'
    version_save = '121821'
    paths = pathsetup('matchingsim')

    filepath = f"{paths['expdatapath']}/{version}/{animal}_all_sessions_{version}.mat"
    fitrangefile = f"{paths['expdatapath']}/102121/fitranges_102121.mat"
    datarange = smart.loadmat(fitrangefile)
    fitrange = datarange['ranges'][datarange['animals'] == animal][0]
    obs, lengths, dirs, fnames, rawchoices = load_multiple_sessions(filepath, fitrange, trialsperblock=15)

    # Find the foraging efficiencies of all blocks
    block_lens = []
    block_corrs = []
    for i in range(len(rawchoices)):
        arr = rawchoices[i]
        block_corrs += list(np.nansum(arr == 1, axis=1))
        block_lens += list(np.sum(~np.isnan(arr), axis=1))

    # Run the fitting procedure
    obs_dim = obs.shape[1]
    np.random.seed(seed)
    masks = ~np.isnan(obs)
    obsmasked = obs[:]
    obsmasked[~masks] = 1

    hmm = ssm.HMM(num_states, obs_dim, observations="blocklapse")
    hmm_lls = hmm.fit(obs, method="em", masks=masks, num_iters=N_iters, init_method="kmeans")

    # Pool states for efficiencies
    zstates = hmm.most_likely_states(obs)
    effs = []
    for i in range(num_states):
        # Find the average foraging efficiency of that state
        blen_state = np.array(block_lens)[zstates == i]
        bcorr_state = np.array(block_corrs)[zstates == i]
        eff_state = sum(bcorr_state) / sum(blen_state)
        effs.append(np.mean(bcorr_state / blen_state))

    # Save the result
    # Save the result
    transmat = hmm.transitions.transition_matrix
    params = hmm.observations.params
    savepath = f"{paths['blockhmmfitpath']}/{version_save}/{animal}_hmmblockfit_{version}_saved_{version_save}.mat"

    vars = ['zstates', 'dirs', 'lengths', 'transmat', 'params',
            'fitrange', 'filepath', 'obs', 'seed', 'hmm_lls', 'effs', 'block_lens', 'block_corrs']
    savedict = make_savedict(vars, locals())

    savefile = 1
    if savefile and not os.path.exists(savepath):
        scipy.io.savemat(savepath, savedict)
        print('File saved')
    elif os.path.exists(savepath):
        print('File exists, skipping save..')


def run_animal(animal, seeds, N_iters=3000, num_states=6):
    '''
    Given animal name and seed number, run block HMM over all the seeds and report the results
    :param animal: str: animal name
    :param seeds: list[int], seed names
    :return: None
    '''
    paths = pathsetup('matchingsim')
    # Load data
    version = '113021'
    filepath = f"{paths['expdatapath']}/{version}/{animal}_all_sessions_{version}.mat"
    fitrangefile = f"{paths['expdatapath']}/102121/fitranges_102121.mat"
    datarange = smart.loadmat(fitrangefile)
    fitrange = datarange['ranges'][datarange['animals'] == animal][0]
    obs, lengths, dirs, fnames, rawchoices = load_multiple_sessions(filepath, fitrange, trialsperblock=15)

    # Run the fitting procedure
    obs_dim = obs.shape[1]
    lls_all = []
    for seed in seeds:
        np.random.seed(seed)
        masks = ~np.isnan(obs)
        obsmasked = obs[:]
        obsmasked[~masks] = 1

        hmm = ssm.HMM(num_states, obs_dim, observations="blocklapse")
        hmm_lls = hmm.fit(obs, method="em", masks=masks, num_iters=N_iters, init_method="kmeans")

        lls_all.append(hmm_lls[-1])
        print(f'animal {animal}, seed value = {seed}, hmm LLS = {hmm_lls[-1]:.2f}')

    # Determine the best seed
    idbest = np2.argmax(lls_all)
    print(f'Best seed is: {seeds[idbest]}')
    return seeds[idbest]

if __name__ == '__main__':
    seeds = [121, 122, 123, 124, 125]
    # animals = ['e46']
    # animals = ['f02', 'f03', 'f04', 'f11', 'f12', 'e35', 'e40',
    #     'fh01', 'fh02', 'f05', 'e53', 'fh03', 'f16', 'f17', 'f20', 'f21', 'f22', 'f23']
    # animals = ['e53', 'fh03', 'f16', 'f17', 'f20', 'f21', 'f22', 'f23']
    paths = pathsetup('matchingsim')
    files = glob.glob(f"{paths['expdatapath']}/113021/*_all_sessions_113021.mat")
    num_states = 6
    N_iters = 3000
    # animals = ['f01']
    animals = [file.split('_')[1].split('/')[-1] for file in files]
    for animal in animals:
        print(f"Running animal: {animal}")
        try:
            bestseed = run_animal(animal, seeds, N_iters, num_states)

            # Run and save with the best seed
            run_and_save(animal, bestseed, N_iters, num_states)
        except:
            continue
