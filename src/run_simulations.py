# from src.utils import simulate_rew_error_correlations, make_switching_world #fit_expfun2, make_switching_world, find_LR_transition_fit, simulate_rew_error_correlations
import src.utils
from src.worldModels import *
from src.agents import *
import numpy as np
import ssm
import tqdm




def run_repeated_single_agent(params):
    '''
    params: dictionary of parameters
    '''
    np.random.seed(params['seed'])

    N_iters = params['N_iters']
    verbose = 1 - params['verbose']

    T11lst = np.zeros(N_iters)
    T22lst = np.zeros(N_iters)
    E1lst = np.zeros(N_iters)
    E2lst = np.zeros(N_iters)
    efflist = np.zeros(N_iters)
    PRslopelist = np.zeros(N_iters)
    PLslopelist = np.zeros(N_iters)
    PRoffsetlist = np.zeros(N_iters)
    PLoffsetlist = np.zeros(N_iters)
    LapseL = np.zeros(N_iters)
    LapseR = np.zeros(N_iters)
    ParamsA = np.zeros(N_iters)
    ParamsB = np.zeros(N_iters)
    ParamsC = np.zeros(N_iters)


    for idx in tqdm.tqdm(range(N_iters), disable=verbose):
        # if idx % 10 == 0:
        #     print('Iteration = ', idx)

        agent, world, pR, pL, hmm = run_single_agent(0, 0, params)

        xvals, means, _, _ = simulate_rew_error_correlations(world, agent)
        paramsFit, _ = fit_expfun2([0.5, 4], xvals, np.array(means))

        efflist[idx] = agent.find_efficiency()

        if params['hmm_fit']:
            T11lst[idx] = hmm.transitions.transition_matrix[0][0]
            T22lst[idx] = hmm.transitions.transition_matrix[1][1]
            E1lst[idx] = logistic(hmm.observations.logit_ps)[0]
            E2lst[idx] = logistic(hmm.observations.logit_ps)[1]

        PRslopelist[idx] = pR[0]
        PLslopelist[idx] = pL[0]
        PRoffsetlist[idx] = pR[1]
        PLoffsetlist[idx] = pL[1]
        LapseL[idx] = pL[2]
        LapseR[idx] = pR[2]
        ParamsA[idx] = paramsFit[0]
        ParamsB[idx] = paramsFit[1]
        # ParamsC[idx] = paramsFit[2]

    return efflist, T11lst, T22lst, E1lst, E2lst, PRslopelist, PLslopelist, \
           PRoffsetlist, PLoffsetlist, LapseL, LapseR, ParamsA, ParamsB, ParamsC

def run_multiple_agents(params):
    '''
    params: dictionary of parameters
    '''
    np.random.seed(params['seed'])

    if params['type']== 'inf-based':
        xlst = params['pswitchlst'] #np.linspace(0.01, 0.45, 10)
        ylst = params['prewlst'] #np.linspace(0.55, 0.99, 10)
    elif params['type'] == 'qlearning':
        xlst = params['gammalst']  # np.linspace(0.01, 0.45, 10)
        ylst = params['epslst']  # np.linspace(0.55, 0.99, 10)
    elif params['type'] == 'v-accumulation':
        xlst = params['gammalst']
        ylst = params['betalst']

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
    ParamsA = np.zeros((len(xlst), len(ylst)))
    ParamsB = np.zeros((len(xlst), len(ylst)))
    ParamsC = np.zeros((len(xlst), len(ylst)))

    coefs_all = []
    perf_train_all = []
    perf_test_all = []

    for idx in range(len(xlst)):
        print('* idx = ', idx)
        for idy in range(len(ylst)):
            print('     idy = ', idy)
            agent, world, pR, pL, hmm = run_single_agent(idx, idy, params)

            coefs, perf_train, perf_test = agent.do_history_logistic_fit(Tmax=params['Tmax'])
            coefs_all.append(coefs)
            perf_train_all.append(perf_train)
            perf_test_all.append(perf_test)

            # Previous fitting code
            xvals, means, _, _ = simulate_rew_error_correlations(world, agent)
            paramsFit, _ = fit_expfun2([0.5, 4], xvals, np.array(means))
            efflist[idx][idy] = agent.find_efficiency()

            if params['hmm_fit']:
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
            ParamsA[idx][idy] = paramsFit[0]
            ParamsB[idx][idy] = paramsFit[1]

    # Reshape the logistic regression results
    coefs_arr = np.array(coefs_all)
    coefs_arr = np.reshape(coefs_arr, [len(xlst), len(ylst), -1])
    perf_train_arr = np.reshape(perf_train_all, [len(xlst), -1])
    perf_test_arr = np.reshape(perf_test_all, [len(xlst), -1])

    return efflist, T11lst, T22lst, E1lst, E2lst, PRslopelist, PLslopelist, \
           PRoffsetlist, PLoffsetlist, LapseL, LapseR, ParamsA, ParamsB, ParamsC, \
            coefs_arr, perf_train_arr, perf_test_arr



def run_single_agent(idx, idy, params):
    '''
    For running a single set of parameters
    '''
    # print(params)
    rlow = params['rlow']
    rhigh = params['rhigh']
    nblocks = params['nblocks']
    ntrials_per_block = params['ntrials_per_block']
    obs_dim = params['obs_dim']
    num_states = params['num_states']
    N_iters = params['N_iters']
    window = params['sigmoid_window']
    # hmm_fit = params['hmm_fit']

    world, _ = src.utils.make_switching_world(rlow, rhigh, nblocks, ntrials_per_block[0], ntrials_per_block[1])

    if params['type'] == 'inf-based':
        agent = EGreedyInferenceBasedAgent(prew=params['prewlst'][idy], pswitch=params['pswitchlst'][idx], eps=params['eps'])
    elif params['type'] == 'qlearning':
        agent = EGreedyQLearningAgent(gamma=params['gammalst'][idx], eps=params['epslst'][idy])
    elif params['type'] == 'v-accumulation':
        agent = ValueAccumulationAgent(gamma=params['gammalst'][idx], beta=params['betalst'][idy])
    else:
        raise ValueError('Agent not recognized, must be inf-based, qlearning, or v-accumulation')

    exp = Experiment(agent, world)
    exp.run()

    # Fit HMM to choice sequence
    if params['hmm_fit']:
        data = np.array(agent.choice_history)[:, None]
        ## testing the constrained transitions class
        hmm = ssm.HMM(num_states, obs_dim, observations="bernoulli")
        _ = hmm.fit(data, method="em", num_iters=N_iters, init_method="kmeans", verbose=0)
    else:
        hmm = None

    # Sigmoidal fit for choice transitions
    # pR, pL, _ = find_LR_transition_fit(world, agent, window=window, type='sigmoid')
    # print('sigmoid:', pR, pL)
    pR, pL, _ = src.utils.find_LR_transition_fit(world, agent, window=window, type='doublesigmoid')
    # print('doublesigmoid', pR, pL)
    return agent, world, pR, pL, hmm

if __name__ == '__main__':
    params = dict(N_iters=50, num_states=2, obs_dim=1, nblocks=100,
                  eps=0, hmm_fit=False, sigmoid_window=30,
                  ntrials_per_block=[5, 20], gammalst=[0.02], epslst=[0.1],
                  rlow=0, rhigh=1, type='qlearning')

    # np.randomuseed(123)
    # world1 = ForagingWorld(prew=0.9, psw=0.1, pstruct=[5, 40], nblockmax=100)

    world = ForagingWorld(prew=0.9, psw=0.1, pstruct=[5, 15], nblockmax=1000)
    # agent = ValueAccumulationAgent(gamma=0.1, beta=10)
    agent = EGreedyQLearningAgent(gamma=0.05, eps=0.01)
    exp = Experiment(agent, world)
    exp.run()
    # agent, world, _, _, _ = run_single_agent(0, 0, params)
    Nerrors = src.utils.get_num_errors_leading_block(world, agent)
    Nrews = src.utils.get_num_rewards_trailing_block(world, agent)
    plt.figure()
    plt.plot(Nrews[:-1], Nerrors[1:], '.')
    plt.show()
    # xvals, means, stds, ysplit = src.utils.simulate_rew_error_correlations(world, agent)
