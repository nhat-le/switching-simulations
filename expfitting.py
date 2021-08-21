import numpy as np
from tqdm import tqdm
import scipy.io
from utils import find_experiment_metrics
import matplotlib.pyplot as plt


def fit_animal(animal, datarange=-1):
    # Load the behavior
    filename = 'expdata/' + animal + '_all_sessions.mat'
    data = scipy.io.loadmat(filename)

    delaystartcands = np.where(data['maxdelays'] > 0)[1]

    if len(delaystartcands) == 0:
        delaystart = len(data['maxdelays'][0])
    else:
        delaystart = delaystartcands[0]
    print('Number of sessions before delay:', delaystart)

    # Fitting for all files
    pL_all = []
    pR_all = []
    eff_all = []

    if datarange == -1:
        datarange = np.arange(0, len(data['targets_cell'][0]))
    elif len(datarange) == 1:
        datarange = np.arange(datarange[0], datarange[0] + 1)
    else:
        datarange = np.arange(datarange[0], datarange[1])

    for i in tqdm(datarange):  # tqdm(range(len(data['targets_cell'][0]))):
        #     print(i)
        expdata = {'alltargets': data['targets_cell'][0][i],
                   'allchoices': data['choices_cell'][0][i]}

        if len(expdata['alltargets']) == 0 or len(np.unique(expdata['alltargets'])) == 1:
            pL_all.append([np.nan] * 3)
            pR_all.append([np.nan] * 3)
            eff_all.append(np.nan)
        else:
            pR, pL, choicelst, eff = find_experiment_metrics(expdata, window=20, type='doublesigmoid')

            # skip if too few blocks
            if choicelst.shape[0] <= 5:
                pL_all.append([np.nan] * 3)
                pR_all.append([np.nan] * 3)
                eff_all.append(np.nan)
            else:
                pL_all.append(pL)
                pR_all.append(pR)
                eff_all.append(eff)

    # Process and return
    pL_all = np.array(pL_all)
    pR_all = np.array(pR_all)
    eff_all = np.array(eff_all)
    slopes = np.vstack([pL_all[:, 0], pR_all[:, 0]])
    offsets = np.vstack([pL_all[:, 1], pR_all[:, 1]])
    lapses = np.vstack([pL_all[:, 2], pR_all[:, 2]])

    offsets[abs(slopes) < 0.01] = -20
    minslopes = np.min(-slopes, axis=0)
    minslopes[minslopes > 1000] = np.nan

    meanoffsets = np.nanmean(-offsets, axis=0)
    maxlapses = np.nanmax(lapses, axis=0)

    return minslopes, meanoffsets, maxlapses, pL_all, pR_all, eff_all

def plot_fitvals(fitparams):
    plt.figure(figsize=(8,4))
    plt.subplot(221)
    plt.plot(fitparams[0], '.')
    plt.title('min slopes')

    plt.subplot(222)
    plt.plot(fitparams[1], '.')
    plt.ylim([0, 30])
    plt.title('mean offsets')

    plt.subplot(223)
    plt.plot(fitparams[2], '.')
    plt.title('max lapses')
    # plt.ylim([0, 10])

    plt.subplot(224)
    plt.plot(fitparams[5], '.')
    plt.title('efficiency')
    plt.tight_layout()

