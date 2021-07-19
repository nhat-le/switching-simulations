from utils import *
from worldModels import *
from agents import *
import numpy as np
import ssm


def run_multiple_agents(params):
    '''
    params: dictionary of parameters
    '''
    np.random.seed(params['seed'])

    if params['type']== 'inf-based':
        pswitchlst = params['pswitchlst'] #np.linspace(0.01, 0.45, 10)
        prewlst = params['prewlst'] #np.linspace(0.55, 0.99, 10)
        xlst = pswitchlst
        ylst = prewlst
    elif params['type'] == 'qlearning':
        gammalst = params['gammalst']  # np.linspace(0.01, 0.45, 10)
        epslst = params['epslst']  # np.linspace(0.55, 0.99, 10)
        xlst = gammalst #TODO: CHECK X/Y VALIDITY
        ylst = epslst

    T11lst = np.zeros((len(xlst), len(ylst)))
    T22lst = np.zeros((len(xlst), len(ylst)))
    E1lst = np.zeros((len(xlst), len(ylst)))
    E2lst = np.zeros((len(xlst), len(ylst)))
    efflist = np.zeros((len(xlst), len(ylst)))
    PRslopelist = np.zeros((len(xlst), len(ylst)))
    PLslopelist = np.zeros((len(xlst), len(ylst)))
    PRoffsetlist = np.zeros((len(xlst), len(ylst)))
    PLoffsetlist = np.zeros((len(xlst), len(ylst)))
    LapseL = np.zeros((len(xlst), len(ylst)))
    LapseR = np.zeros((len(xlst), len(ylst)))

    for idx in range(len(xlst)):
        print('* idx = ', idx)
        for idy in range(len(ylst)):
            print('     idy = ', idy)

            agent, world, pR, pL, hmm = run_single_agent(idx, idy, params)

            efflist[idx][idy] = agent.find_efficiency()
            T11lst[idx][idy] = hmm.transitions.transition_matrix[0][0]
            T22lst[idx][idy] = hmm.transitions.transition_matrix[1][1]
            E1lst[idx][idy] = logistic(hmm.observations.logit_ps)[0]
            E2lst[idx][idy] = logistic(hmm.observations.logit_ps)[1]
            PRslopelist[idx][idy] = pR[0]
            PLslopelist[idx][idy] = pL[0]
            PRoffsetlist[idx][idy] = pR[1]
            PLoffsetlist[idx][idy] = pL[1]
            LapseL[idx][idy] = pL[2]
            LapseR[idx][idy] = pR[2]

    return efflist, T11lst, T22lst, E1lst, E2lst, PRslopelist, PLslopelist, \
           PRoffsetlist, PLoffsetlist, LapseL, LapseR



def run_single_agent(idx, idy, params):
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

    if params['type'] == 'inf-based':
        agent = EGreedyInferenceBasedAgent(prew=params['prewlst'][idy], pswitch=params['pswitchlst'][idx], eps=params['eps'])
    elif params['type'] == 'qlearning':
        agent = EGreedyQLearningAgent(gamma=params['gammalst'][idx], eps=params['epslst'][idy])

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

