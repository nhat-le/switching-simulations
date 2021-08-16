from run_simulations import *
from sklearn import svm



def get_Qmetrics(gamma, eps):
    gammalst = [gamma]
    epslst = [eps]

    agent_type = 'qlearning'  # 'qlearning' or 'inf-based' or 'v-accumulation'

    N_iters = 50
    num_states = 2
    obs_dim = 1
    nblocks = 20  # 100
    eps = 0
    hmm_fit = False
    sigmoid_window = 30
    ntrials_per_block = [10, 40]
    seed = 0
    rlow = 0
    rhigh = 1
    np.random.seed(123)

    params = {'N_iters': N_iters, 'num_states': num_states, 'obs_dim': obs_dim,
              'nblocks': nblocks, 'eps': eps, 'ntrials_per_block': ntrials_per_block,
              'gammalst': gammalst, 'epslst': epslst, 'seed': seed, 'type': agent_type, 'hmm_fit': hmm_fit,
              'seed': seed, 'sigmoid_window': sigmoid_window, 'rlow': rlow, 'rhigh': rhigh, 'verbose': False}

    # Q-learning simulation
    #     print('Starting q-learning simulation')
    efflistQ, T11lstQ, T22lstQ, E1lstQ, E2lstQ, PRslopelistQ, PLslopelistQ, \
    PRoffsetlistQ, PLoffsetlistQ, LapseLQ, LapseRQ, ParamsAQ, ParamsBQ, ParamsCQ = run_repeated_single_agent(params)

    Qmetrics = [efflistQ, LapseLQ, PLoffsetlistQ, PLslopelistQ]
    return Qmetrics


def get_IB_metrics(psw, prew):
    pswitchlst = [psw]
    prewlst = [prew]

    agent_type = 'inf-based'  # 'qlearning' or 'inf-based' or 'v-accumulation'

    N_iters = 50
    num_states = 2
    obs_dim = 1
    nblocks = 20  # 100
    eps = 0
    hmm_fit = False
    sigmoid_window = 30
    ntrials_per_block = [10, 40]
    seed = 0
    rlow = 0
    rhigh = 1
    np.random.seed(123)

    params = {'N_iters': N_iters, 'num_states': num_states, 'obs_dim': obs_dim,
              'nblocks': nblocks, 'eps': eps, 'ntrials_per_block': ntrials_per_block,
              'seed': seed, 'type': agent_type, 'hmm_fit': hmm_fit,
              'pswitchlst': pswitchlst, 'prewlst': prewlst,
              'seed': seed, 'sigmoid_window': sigmoid_window, 'rlow': rlow, 'rhigh': rhigh, 'verbose': False}

    # Inference-based simulation
    efflistIB, T11lstIB, T22lstIB, E1lstIB, E2lstIB, PRslopelistIB, PLslopelistIB, \
    PRoffsetlistIB, PLoffsetlistIB, LapseLIB, LapseRIB, ParamsAIB, ParamsBIB, ParamsCIB = run_repeated_single_agent(
        params)

    IBmetrics = [efflistIB, LapseLIB, PLoffsetlistIB, PLslopelistIB]
    return IBmetrics


def fit_and_evaluate_svm(Qmetrics, IBmetrics):
    # Make the training and test sets
    XQ = np.vstack(Qmetrics)
    XIB = np.vstack(IBmetrics)
    Xall = np.hstack([XQ, XIB])
    yall = np.array([0] * XQ.shape[1] + [1] * XIB.shape[1])

    # Normalize
    mean = np.mean(Xall, axis=1)
    std = np.std(Xall, axis=1)
    Xall = (Xall.T - mean) / std

    # Shuffle
    order = np.arange(Xall.shape[0])
    np.random.shuffle(order)
    Xall = Xall[order, :]

    Xtrain = Xall[:75, :]
    Xtest = Xall[75:, :]
    ytrain = yall[order][:75]
    ytest = yall[order][75:]

    clf = svm.SVC(kernel='linear')
    clf.fit(Xtrain, ytrain)
    coefs = clf.coef_

    preds = clf.predict(Xtest)
    perf = sum(preds == ytest) / len(preds)
    return perf, coefs


def find_svm_perf(gamma, eps, psw, prew):
    pswitchlst = [psw]
    prewlst = [prew]

    gammalst = [gamma]
    epslst = [eps]

    agent_type = 'qlearning'  # 'qlearning' or 'inf-based' or 'v-accumulation'

    N_iters = 50
    num_states = 2
    obs_dim = 1
    nblocks = 20  # 100
    eps = 0
    hmm_fit = False
    sigmoid_window = 30
    ntrials_per_block = [10, 40]
    seed = 0
    rlow = 0
    rhigh = 1
    np.random.seed(123)

    params = {'N_iters': N_iters, 'num_states': num_states, 'obs_dim': obs_dim,
              'nblocks': nblocks, 'eps': eps, 'ntrials_per_block': ntrials_per_block,
              'gammalst': gammalst, 'epslst': epslst, 'seed': seed, 'type': agent_type, 'hmm_fit': hmm_fit,
              'pswitchlst': pswitchlst, 'prewlst': prewlst,
              'seed': seed, 'sigmoid_window': sigmoid_window, 'rlow': rlow, 'rhigh': rhigh, 'verbose' :False}

    # Q-learning simulation
    print('Starting q-learning simulation')
    efflistQ, T11lstQ, T22lstQ, E1lstQ, E2lstQ, PRslopelistQ, PLslopelistQ, \
    PRoffsetlistQ, PLoffsetlistQ, LapseLQ, LapseRQ, ParamsAQ, ParamsBQ, ParamsCQ = run_repeated_single_agent(params)

    # Inference-based simulation
    print('Starting inf-based simulation')
    params['type'] = 'inf-based'
    efflistIB, T11lstIB, T22lstIB, E1lstIB, E2lstIB, PRslopelistIB, PLslopelistIB, \
    PRoffsetlistIB, PLoffsetlistIB, LapseLIB, LapseRIB, ParamsAIB, ParamsBIB, ParamsCIB = run_repeated_single_agent \
        (params)

    Qmetrics = [efflistQ, LapseLQ, PLoffsetlistQ, PLslopelistQ]
    IBmetrics = [efflistIB, LapseLIB, PLoffsetlistIB, PLslopelistIB]
    perf, coefs = fit_and_evaluate_svm(Qmetrics, IBmetrics)
    return Qmetrics, IBmetrics, perf, coefs


def normalize_and_average(arr):
    normarr = (arr - np.mean(arr)) / np.std(arr)
    return np.mean(normarr, axis=2)

def normalize_and_average_all(lst):
    return [normalize_and_average(elem) for elem in lst]


