import numpy as np
from worldModels import *
# from agents import *
import scipy.optimize

# import autograd.numpy as np
# import autograd.numpy.random as npr
# npr.seed(0)

import ssm
from ssm.util import find_permutation
from ssm.plots import gradient_cmap, white_to_color_cmap

def logistic(x, beta=1):
    return 1 / (1 + np.exp(-beta * x))


def make_switching_world(rlow, rhigh, nblocks, ntrialsLow, ntrialsHigh):
    ratesL = (np.mod(np.arange(nblocks), 2)) * (rhigh - rlow) + rlow
    ratesR = (1 - np.mod(np.arange(nblocks), 2)) * (rhigh - rlow) + rlow
    if np.random.rand() > 0.5:
        rates = np.vstack((ratesL, ratesR)).T
    else:
        rates = np.vstack((ratesR, ratesL)).T
    ntrials = np.random.uniform(low=ntrialsLow, high=ntrialsHigh, size=nblocks).astype('int')
    world = PersistentWorld(rates=rates, ntrials=ntrials)
    return world, ntrials


def make_switching_world_withCheck(rlow, rhigh, nblocks, ntrialsLow, ntrialsHigh):
    ratesL = (np.mod(np.arange(nblocks), 2)) * (rhigh - rlow) + rlow
    ratesR = (1 - np.mod(np.arange(nblocks), 2)) * (rhigh - rlow) + rlow
    if np.random.rand() > 0.5:
        rates = np.vstack((ratesL, ratesR)).T
    else:
        rates = np.vstack((ratesR, ratesL)).T
    ntrials = np.random.uniform(low=ntrialsLow, high=ntrialsHigh, size=nblocks).astype('int')
    world = PersistentWorldWithCheck(rates=rates, ntrials=ntrials, threshold=0.8)
    return world, ntrials


def errorsigmoid(p, x, y):
    '''
    Error function used for sigmoid fitting
    '''
    lapse = p[2]
    preds = lapse + (1 - 2 * lapse) * 1 / (1 + np.exp(-p[0] * (x + p[1])))

    return np.sum((preds - y) ** 2)


def find_LR_transition_fit(world, agent, window):
    '''
    For a switching world, determines the agent transition functions,
    for left->right and right->left transitions
    window: how many trials after the transition do we want to keep for fitting?
    '''
    # Find where the block transitions happen
    # side_history = np.array(world.rate_history)[:, 0]
    # blocktrans = np.where(np.diff(side_history) != 0)[0][:-1]

    # Get the choices around the transition
    if world.ntrialblocks[-1] == 0:
        choicelst = split_by_trials(np.array(agent.choice_history), world.ntrialblocks[:-1], chop='min')
    else:
        choicelst = split_by_trials(np.array(agent.choice_history), world.ntrialblocks, chop='min')

    if window is None:
        window = choicelst.shape[1]
    choicelst = choicelst[:, :window]
    # print(choicelst.shape)


    #     choicelst = []
    #     for i in range(window):
    #         choicelst.append(np.array(agent.choice_history)[blocktrans + i])

    #     choicelst = np.array(choicelst)

    # print('choicemean = ', np.mean(choicelst[:,::2]), 'side=  ', world.side_history[0][0])
    if np.ndim(np.array(world.side_history)) == 1:
        # print('here, first side = ', )
        pRight, pLeft = fit_sigmoidal(choicelst, first_side=world.side_history[0])
    else:
        # print('here 2, first side = ', world.side_history[0][0])
        # pRight, pLeft = fit_sigmoidal(choicelst, first_side=world.side_history[0][0])
        pRight, pLeft = fit_sigmoidal(choicelst, first_side=world.rate_history[0][0] < 0.5)
    return (pRight, pLeft)


def find_transition_guess(sig):
    '''
    Returns the offset where a time series crosses 0.5
    '''
    return np.argmin((sig - 0.5) ** 2)


def find_transition_guess_binary(sig):
    '''
    Returns the offset where a time series crosses 0.5, through binary segmentation
    '''
    candidates = np.where(np.diff(sig > 0.5) != 0)[0]
    if len(candidates) == 0:
        return -1
    else:
        return np.where(np.diff(sig > 0.5) != 0)[0][0]


def fit_sigmoidal(choicelst, first_side):
    '''
    Fit a sigmoidal to the average choice data
    first_side: first side that is rewarded, i.e. world.side_history[0][0]
    returns: pright, pleft, where each is a tuple (slope, offset, lapse)
    '''
    # print('choicemean = ', np.mean(choicelst[:,::2]), 'side=  ', first_side)
    if first_side == 0:
        # print('left')
        leftAverage = np.nanmean(choicelst[1::2, :], axis=0)
        rightAverage = np.nanmean(choicelst[::2, :], axis=0)
    else:
        # print('right')
        rightAverage = np.nanmean(choicelst[1::2, :], axis=0)
        leftAverage = np.nanmean(choicelst[::2, :], axis=0)

    offsetsR = np.arange(len(rightAverage))
    offsetsL = np.arange(len(leftAverage))

    # Fit right transitions
    # print(rightAverage)
    funR = lambda x: errorsigmoid(x, offsetsR, rightAverage)
    switchGuessR = find_transition_guess_binary(rightAverage)  # offset that crosses 0.5
    if switchGuessR == -1:  # No switch happened!
        # pRight = [0, -np.inf, 0]
        paramsRight = scipy.optimize.minimize(funR, [1, -len(rightAverage), 0],
                                              bounds=((None, 0), (None, 0), (0, 0.5)))
    else:
        paramsRight = scipy.optimize.minimize(funR, [1, -switchGuessR, 0],
                                              bounds=((None, 0), (None, 0), (0, 0.5)))
    pRight = paramsRight.x
    # print(pRight)

    # print('done with right')
    # Fit left transitions
    # print(leftAverage)
    funL = lambda x: errorsigmoid(x, offsetsL, leftAverage)
    switchGuessL = find_transition_guess_binary(leftAverage)
    if switchGuessL == -1:  # No switch happened!
        # pLeft = [0, -np.inf, 0]
        paramsLeft = scipy.optimize.minimize(funL, [-1, -len(leftAverage), 0],
                                             bounds=((0, None), (None, 0), (0, 0.5)))
    else:
        # print('here')
        paramsLeft = scipy.optimize.minimize(funL, [-1, -switchGuessL, 0],
                                             bounds=((0, None), (None, 0), (0, 0.5)))
    pLeft = paramsLeft.x
    # print(pLeft)
    # print('done with left')

    return pRight, pLeft


def split_by_trials(seq, ntrials, chop='none'):
    '''
    seq: an array of integers, of length sum(ntrials)
    ntrials: number of trials per block
    splits the seq into small chunks with lengths specified by ntrials
    chop: if none, no chopping, if min, chop to the shortest length (min(ntrials)),
    if max, pad to the longest length (min(ntrials)),
    '''
    if ntrials[-1] == 0:
        ntrials = ntrials[:-1]

    minN = min(ntrials)
    maxN = max(ntrials)
    if len(seq) != sum(ntrials):
        raise ValueError('ntrials must sum up to length of sequence')
    endpoints = np.cumsum(ntrials)[:-1]

    splits = np.split(seq, endpoints)

    if chop == 'min':
        # print('here')
        # print(np.array([elem[:minN] for elem in splits]))
        return np.array([elem[:minN] for elem in splits])

    elif chop == 'none':
        return splits
    elif chop == 'max':
        # pad to the max len
        result = np.ones((len(ntrials), maxN)) * np.nan
        for i in range(len(ntrials)):
            result[i, 0:ntrials[i]] = splits[i]
        return result
    else:
        raise ValueError('invalide chop type')


def get_zstates(agent, num_states=2, obs_dim=1, N_iters=50):
    '''
    Fit HMM to the choice sequence of the agent
    Returns: the sequence of the most likely z-states
    '''
    # Fit HMM to choice sequence
    data = np.array(agent.choice_history)[:, None]

    ## testing the constrained transitions class
    hmm = ssm.HMM(num_states, obs_dim, observations="bernoulli")
    hmm_lls = hmm.fit(data, method="em", num_iters=N_iters, init_method="kmeans", verbose=0)

    if hmm.observations.logit_ps[0] > 0:
        return 1 - hmm.most_likely_states(data)
    else:
        return hmm.most_likely_states(data)


def get_switch_times(world, agent):
    '''
    Returns an array of switching times (in trials),
    based on the HMM model fits
    '''
    z_states = get_zstates(agent)
    splits = split_by_trials(z_states, world.ntrialblocks, chop='none')
    # Identify where the switch happens
    if np.ndim(np.array(world.side_history)) == 1:
        first_side = world.side_history[0]
    else:
        first_side = world.side_history[0][1]

    switchlst = []
    for i in range(len(splits)):
        arr = splits[i]
        print(arr)
        print(i, first_side)
        # Skip trials that start on the wrong side
        if arr[0] == (first_side + i) % 2:
            switch = -1
            print('skipping')
        else:
            # Find the first element that is the opposite state
            target = (i + first_side) % 2
            cands = np.where(arr == target)[0]
            if len(cands) == 0:
                switch = world.ntrialblocks[i]
            else:
                switch = cands[0]
        # print('switch =', switch)
        switchlst.append(switch)

    return np.array(switchlst)


def get_num_errors_leading_block(world, agent):
    '''
    Returns an array of switching times (in trials),
    switch times based on the first time animal switches
    '''
    choicelst = split_by_trials(agent.choice_history, world.ntrialblocks, chop='none')

    if np.ndim(np.array(world.side_history)) == 1:
        first_side = world.side_history[0]
    else:
        first_side = world.side_history[0][0]

    switchlst = []
    for i in range(len(world.ntrialblocks) - 1):
        blockchoice = choicelst[i]
        target = (first_side + i) % 2
        if blockchoice[0] == target:
            switch = -1
        else:
            switch = np.where(blockchoice == target)[0][0]

        switchlst.append(switch)

    return np.array(switchlst)


def get_num_rewards_trailing_block(world, agent):
    '''
    Returns an array of consec rewards at the end of the block
    switch times based on the first time animal switches
    '''
    choicelst = split_by_trials(agent.choice_history, world.ntrialblocks, chop='none')

    if np.ndim(np.array(world.side_history)) == 1:
        first_side = world.side_history[0]
    else:
        first_side = world.side_history[0][0]

    nrewlst = []
    for i in range(len(world.ntrialblocks) - 1):
        blockchoice = choicelst[i]
        blockchoiceflip = np.flip(blockchoice)
        target = (first_side + i) % 2
        #         print('array is ', blockchoice)
        if blockchoiceflip[0] != target:
            nrew = -1
        #             print('skipping')
        else:
            nrew = np.where(blockchoiceflip != target)[0]
            #             print(nrew)
            if len(nrew) == 0:
                nrew = len(blockchoiceflip)
            else:
                nrew = nrew[0]
        #             print(target)
        #             print(nrew)

        nrewlst.append(nrew)

    return np.array(nrewlst)


def exp_fun(x, params):
    alpha = params[0]
    beta = params[1]
    C = params[2]
    return C - beta * np.exp(-alpha * x)


def loss(params, x, y):
    pred = exp_fun(x, params)
    return np.sum((pred - y) ** 2)


def fit_expfun(params0, datax, datay):
    # Filter out nan's in datax and datay
    goody = datay[~np.isnan(datax)]
    goodx = datax[~np.isnan(datax)]

    result = scipy.optimize.minimize(loss, params0, (goodx, goody))
    params = result.x
    ypred = exp_fun(datax, params)
    return params, ypred

def simulate_rew_error_correlations(world, agent):
    exp = Experiment(agent, world)
    exp.run()

    lst = get_switch_times(world, agent).astype('float')
    print(lst)
    lst[lst == -1] = np.nan
    nafterswitch = world.ntrialblocks[:-1] - lst

    # Aggregate data
    xarr = nafterswitch[:-1]
    yarr = lst[1:]
    order = np.argsort(xarr)
    xsorted = xarr[order]
    ysortbyX = yarr[order]

    xvals, idx = np.unique(xsorted, return_index=True)
    # print(xsorted)
    ysplit = np.split(ysortbyX, idx[1:])
    # print(ysplit)

    # Mean of each split
    means = []
    stds = []
    for elem in ysplit[1:]:
        # print(elem)
        means.append(np.nanmean(elem))
        stds.append(np.nanstd(elem) / np.sqrt(len(elem)))

    return xvals[1:], means, stds, ysplit




