from utils import *
from worldModels import *
import numpy as np
import ssm
from ssm.util import find_permutation
from ssm.plots import gradient_cmap, white_to_color_cmap


def run_inf_based(params):
    '''
    params: dictionary of parameters
    '''
    np.random.seed(params['seed'])
    pswitchlst = params['pswitchlst'] #np.linspace(0.01, 0.45, 10)
    prewlst = params['prewlst'] #np.linspace(0.55, 0.99, 10)
    eps = params['eps']

    T11lst = np.zeros((len(pswitchlst), len(prewlst)))
    T22lst = np.zeros((len(pswitchlst), len(prewlst)))
    E1lst = np.zeros((len(pswitchlst), len(prewlst)))
    E2lst = np.zeros((len(pswitchlst), len(prewlst)))
    efflist = np.zeros((len(pswitchlst), len(prewlst)))
    PRslopelist = np.zeros((len(pswitchlst), len(prewlst)))
    PLslopelist = np.zeros((len(pswitchlst), len(prewlst)))
    PRoffsetlist = np.zeros((len(pswitchlst), len(prewlst)))
    PLoffsetlist = np.zeros((len(pswitchlst), len(prewlst)))
    LapseL = np.zeros((len(pswitchlst), len(prewlst)))
    LapseR = np.zeros((len(pswitchlst), len(prewlst)))

    for idsw, psw in enumerate(pswitchlst):
        print('* pswitch = ', psw)
        for idrew, prew in enumerate(prewlst):
            print('     prew = ', prew)
            agent, world, pR, pL, hmm = run_single_inf_based_agent(prew, psw, eps, params)

            efflist[idsw][idrew] = agent.find_efficiency()
            T11lst[idsw][idrew] = hmm.transitions.transition_matrix[0][0]
            T22lst[idsw][idrew] = hmm.transitions.transition_matrix[1][1]
            E1lst[idsw][idrew] = logistic(hmm.observations.logit_ps)[0]
            E2lst[idsw][idrew] = logistic(hmm.observations.logit_ps)[1]
            PRslopelist[idsw][idrew] = pR[0]
            PLslopelist[idsw][idrew] = pL[0]
            PRoffsetlist[idsw][idrew] = pR[1]
            PLoffsetlist[idsw][idrew] = pL[1]
            LapseL[idsw][idrew] = pL[2]
            LapseR[idsw][idrew] = pR[2]

    return efflist, T11lst, T22lst, E1lst, PRslopelist, PLslopelist, \
           PRoffsetlist, PLoffsetlist, LapseL, LapseR



def run_single_inf_based_agent(prew, psw, eps, params):
    '''
    For running a single set of parameters
    '''
    rlow = params['rlow']
    rhigh = params['rhigh']
    nblocks = params['nblocks']
    ntrials_per_block = params['ntrials_per_block']
    obs_dim = params['obs_dim']
    num_states = params['num_states']
    N_iters = params['N_iters']

    world, _ = make_switching_world(rlow, rhigh, nblocks, ntrials_per_block[0], ntrials_per_block[1])
    agent = EGreedyInferenceBasedAgent(prew=prew, pswitch=psw, eps=eps)
    exp = Experiment(agent, world)
    exp.run()

    # Fit HMM to choice sequence
    data = np.array(agent.choice_history)[:, None]

    ## testing the constrained transitions class
    hmm = ssm.HMM(num_states, obs_dim, observations="bernoulli")
    _ = hmm.fit(data, method="em", num_iters=N_iters, init_method="kmeans", verbose=0)

    # Sigmoidal fit for choice transitions
    pR, pL = find_LR_transition_fit(world, agent, window=15)
    return agent, world, pR, pL, hmm
