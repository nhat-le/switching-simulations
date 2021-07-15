import numpy as np
from worldModels import *
import scipy.optimize

def logistic(x):
    return 1 / (1 + np.exp(-x))


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
    side_history = np.array(world.rate_history)[:, 0]
    blocktrans = np.where(np.diff(side_history) != 0)[0][:-1]

    # Get the choices around the transition
    choicelst = split_by_trials(np.array(agent.choice_history), world.ntrialblocks, chop='max')
    #     choicelst = []
    #     for i in range(window):
    #         choicelst.append(np.array(agent.choice_history)[blocktrans + i])

    #     choicelst = np.array(choicelst)

    # print('choicemean = ', np.mean(choicelst[:,::2]), 'side=  ', world.side_history[0][0])
    if np.ndim(np.array(world.side_history)) == 1:
        pRight, pLeft = fit_sigmoidal(choicelst, first_side=world.side_history[0])
    else:
        pRight, pLeft = fit_sigmoidal(choicelst, first_side=world.side_history[0][0])
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
    '''
    # print('choicemean = ', np.mean(choicelst[:,::2]), 'side=  ', first_side)
    if first_side == 0:
        # print('left')
        leftAverage = np.mean(choicelst[:, ::2], axis=1)
        rightAverage = np.mean(choicelst[:, 1::2], axis=1)
    else:
        # print('right')
        rightAverage = np.mean(choicelst[:, ::2], axis=1)
        leftAverage = np.mean(choicelst[:, 1::2], axis=1)

    offsets = np.arange(len(leftAverage))

    # Fit right transitions
    funR = lambda x: errorsigmoid(x, offsets, rightAverage)
    switchGuessR = find_transition_guess_binary(rightAverage)  # offset that crosses 0.5
    if switchGuessR == -1:  # No switch happened!
        pRight = [0, -np.inf, 0]
    else:
        paramsRight = scipy.optimize.minimize(funR, [1, -switchGuessR, 0])
        pRight = paramsRight.x

    funL = lambda x: errorsigmoid(x, offsets, leftAverage)
    switchGuessL = find_transition_guess_binary(leftAverage)
    if switchGuessL == -1:  # No switch happened!
        pLeft = [0, -np.inf, 0]
    else:
        paramsLeft = scipy.optimize.minimize(funL, [-1, -switchGuessL, 0])
        pLeft = paramsLeft.x

    return pRight, pLeft


def split_by_trials(seq, ntrials, chop='none'):
    '''
    seq: an array of integers, of length sum(ntrials)
    ntrials: number of trials per block
    splits the seq into small chunks with lengths specified by ntrials
    chop: if none, no chopping, if min, chop to the shortest length (min(ntrials)),
    if max, pad to the longest length (min(ntrials)),
    '''
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