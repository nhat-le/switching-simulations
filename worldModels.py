import numpy as np

class Agent():
    '''
    A general agent
    '''
    def __init__(self):
        self.outcome_history = []
        self.choice_history = []
        self.Rewards1side_history = []
        self.Rewards0side_history = []
        
    def find_efficiency(self):
        return sum(self.outcome_history) / len(self.outcome_history)

    def find_prob(self):
        '''
        Find the instantaneous probability of the agent for all trials
        '''
        p0 = np.array(self.Rewards0side_history)
        p1 = np.array(self.Rewards1side_history)
        return p1 / (p0 + p1)

    def get_running_choice_fraction(self, window):
        '''
        Get the running mean of choice fraction
        window: int, window for averaging
        '''
        choice_arr = np.array(self.choice_history)
        kernel = np.ones(window) / window
        return np.convolve(choice_arr, kernel, 'same')


    def get_running_reward_fraction(self, window):
        '''
        Get the running mean of reward fraction
        window: int, window for averaging
        '''
        rewards0side = np.array(self.Rewards0side_history)
        rewards1side = np.array(self.Rewards1side_history)
        kernel = np.ones(window) / window

        running_0sidefrac = np.convolve(rewards0side, kernel, 'same')
        running_1sidefrac = np.convolve(rewards1side, kernel, 'same')

        return running_1sidefrac / (running_0sidefrac + running_1sidefrac)

class World():
    '''
    A general class for worlds
    '''
    def __init__(self):
        self.history = []
        self.rate_history = []

    def find_prob(self):
        '''
        Returns array of rates for trial-1 side
        '''
        rateArr = np.array(self.rate_history)
        return rateArr[:,1]

    def get_rate_history(self):
        return np.array(self.rate_history)


# A world object
class RandomWorld(World):
    '''
    A probablistic world that alternates between blocks of constant probability
    '''
    def __init__(self, rates, ntrials):
        '''
        rates: a list, rates of the world
        ntrials: a list, number of trials in each block
        '''
        self.rates = rates
        self.ntrials = ntrials
        self.curr_block = 0
        self.history = []
        self.rate_history = []
        self.curr_rate = self.rates[0]
        self.curr_side = np.random.rand() < self.curr_rate
        #print('curr side=  ', int(self.curr_side))
        
    def update(self, agent_choice):
        '''
        Update the world based on agent choice
        '''
        self.history.append(self.curr_side)
        self.rate_history.append(self.curr_rate)
        
        if agent_choice == self.curr_side:
            reward = 1
            #print('len history = ', len(self.history), 'next switch at', sum(self.ntrials[:self.curr_block + 1]))
            # See if the world should switch blocks
            if len(self.history) > sum(self.ntrials[:self.curr_block + 1]):
#                 print('world switching!')
                self.curr_block += 1
                self.curr_rate = self.rates[self.curr_block]
                        
            # Sample the next side
            self.curr_side = np.random.rand() < self.curr_rate
        
        else:
            # Incorrect, keep the same side
            reward = 0
                    
        return reward
    

# A world object
class PersistentWorld(World):
    '''
    A probablistic world that alternates between blocks of constant probability
    '''
    def __init__(self, rates, ntrials):
        '''
        rates: a nblocks x 2 array, representing rates for 0- and 1-sides
        ntrials: a list, number of trials in each block
        '''
        self.rates = rates
        self.ntrials = ntrials
        self.curr_block = 0
        self.side_history = []
        self.rate_history = []
        
        self.curr_rates = np.array(self.rates[0,:])        
        self.active_sites = np.random.rand(2) < rates[0,:]
        
#         print('curr active sites:', self.active_sites)
        
        #print('curr side=  ', int(self.curr_side))
        
    def update(self, agent_choice):
        '''
        Update the world based on agent choice
        '''
        agent_choice = int(agent_choice)
        self.rate_history.append(self.curr_rates.copy())
        self.side_history.append(self.active_sites.copy())
        
        # Is there reward at current choice side?
        reward = self.active_sites[agent_choice]
#         print('choice = ', agent_choice, 'reward = ', reward)
#         print('n trials so far =', len(self.side_history))
        
        
        # Are we switching blocks?
        if len(self.rate_history) > sum(self.ntrials[:self.curr_block + 1]):
            self.curr_block += 1
            self.curr_rates = self.rates[self.curr_block,:]
#             print('world switching! curr rates = ', self.curr_rates, 'trials so far =', len(self.side_history))
        
        
        # Update active_sites
        if self.active_sites[0] == 0 or agent_choice == 0:
#             print('updated site 0')
            self.active_sites[0] = np.random.rand() < self.curr_rates[0]

        if self.active_sites[1] == 0 or agent_choice == 1:
#             print('updated site 1', self.curr_rates[1])
            self.active_sites[1] = np.random.rand() < self.curr_rates[1]

#         print('current active sites: ', self.active_sites)
            
                    
        return reward


class PersistentWorldWithCheck(World):
    '''
    A probablistic world that alternates between blocks of constant probability,
    only switching after performance crosses a threshold
    '''

    def __init__(self, rates, ntrials, threshold):
        '''
        rates: a nblocks x 2 array, representing rates for 0- and 1-sides
        ntrials: a list, number of trials in each block
        '''
        self.rates = rates
        self.threshold = threshold
        self.currperf = -1 # Current performance in the block

        self.currCorrect = 0
        self.currIncorr = 0

        self.ntrials = ntrials
        self.curr_block = 0
        self.side_history = []
        self.rate_history = []

        self.curr_rates = np.array(self.rates[0, :])
        self.active_sites = np.random.rand(2) < rates[0, :]

    def update(self, agent_choice):
        '''
        Update the world based on agent choice
        '''
        agent_choice = int(agent_choice)
        self.rate_history.append(self.curr_rates.copy())
        self.side_history.append(self.active_sites.copy())

        # Is there reward at current choice side?
        reward = self.active_sites[agent_choice]

        if reward > 0:
            self.currCorrect += 1
        else:
            self.currIncorr += 1

        self.currperf = self.currCorrect / (self.currCorrect + self.currIncorr)
        #print(self.currperf)


        # Are we switching blocks?
        if self.currCorrect + self.currIncorr > self.ntrials[self.curr_block] and \
            self.currperf > self.threshold:
            #print('Block switching!')
            self.curr_block += 1
            self.curr_rates = self.rates[self.curr_block, :]
            self.currCorrect = 0
            self.currIncorr = 0
            self.currperf = -1

        # Update active_sites
        if self.active_sites[0] == 0 or agent_choice == 0:
            #             print('updated site 0')
            self.active_sites[0] = np.random.rand() < self.curr_rates[0]

        if self.active_sites[1] == 0 or agent_choice == 1:
            #             print('updated site 1', self.curr_rates[1])
            self.active_sites[1] = np.random.rand() < self.curr_rates[1]

        #         print('current active sites: ', self.active_sites)

        return reward

            
    
class MatchingAgent(Agent):
    '''
    Simulate an agent that matches perfectly (perfect integration of past rewards)
    '''
    def __init__(self, eps):
        '''
        eps: rate is limited to the range [eps, 1-eps]
        '''
        assert(0 <= eps <= 0.5)
        self.Rewards1side = 1 # Start with 1 so that the agent does not get 'stuck'
        self.Rewards0side = 1
        self.eps = eps
        self.choice_history = []
        self.Rewards1side_history = []
        self.Rewards0side_history = []
        self.outcome_history = []
        
    def outcome_received(self, outcome):
        if outcome == 1:
            if self.choice_history[-1] == 1:
                self.Rewards1side += 1
            else:
                self.Rewards0side += 1
                
        self.outcome_history.append(outcome)
        self.Rewards1side_history.append(self.Rewards1side)
        self.Rewards0side_history.append(self.Rewards0side)
              
    
    def make_choice(self):
        '''
        Make a choice, probabilistically sample from past reward ratios
        '''
        
        p = self.Rewards1side / (self.Rewards1side + self.Rewards0side)
        p = min(p, 1 - self.eps)
        p = max(p, self.eps)
        choice = np.random.rand() < p
        self.choice_history.append(choice)
        return choice


class EGreedyQLearningAgent(Agent):
    '''
    Simulate an agent that uses Q-learning, and chooses action probabilistically
    based on the ratio of the q-values
    '''

    def __init__(self, gamma, eps=0):
        '''
        gamma: learning rate for q-value updates
        eps: used for epsilon-greedy strategy, eps = 0 means a greedy agent
        '''
        self.q0 = 0.5
        self.q1 = 0.5
        self.gamma = gamma
        self.eps = eps
        self.choice_history = []
        self.q0_history = []
        self.q1_history = []
        self.outcome_history = []
        self.Rewards0side_history = []
        self.Rewards1side_history = []

    def outcome_received(self, outcome):
        if self.choice_history[-1] == 1:
            self.q1 = self.q1 + self.gamma * (outcome - self.q1)
        else:
            self.q0 = self.q0 + self.gamma * (outcome - self.q0)

        self.outcome_history.append(outcome)
        self.q0_history.append(self.q0)
        self.q1_history.append(self.q1)

        # Update reward history on each side
        if outcome == 1:
            if self.choice_history[-1] == 1:
                self.Rewards1side_history.append(1)
                self.Rewards0side_history.append(0)
            else:
                self.Rewards1side_history.append(0)
                self.Rewards0side_history.append(1)

    def make_choice(self):
        '''
        Make a choice using the epsilon-greedy strategy
        '''
        # Flip a coin to decide if explore or exploit
        explore = np.random.rand() < self.eps
        if explore: #choose actions randomly
            choice = np.random.rand() < 0.5
        else:
            if self.q1 > self.q0:
                choice = 1
            elif self.q1 < self.q0:
                choice = 0
            else:
                choice = np.random.rand() < 0.5

        self.choice_history.append(choice)
        return choice

    def find_prob(self):
        '''
        Find the instantaneous probability of the agent for all trials
        '''
        p0 = np.array(self.q0_history)
        p1 = np.array(self.q1_history)
        return p1 / (p0 + p1)
        
    
# A Matching agent object
class ConstantProbAgent(Agent):
    '''
    Simulate an agent that decides with a fixed probability
    '''
    def __init__(self, prob):
        self.choice_history = []
        self.prob = prob
        self.outcome_history = []
#         self.Rewards1side_history = []
#         self.Rewards0side_history = []
        
    def outcome_received(self, outcome):
        self.outcome_history.append(outcome)
        
    
    def make_choice(self):
        '''
        Make a choice, probabilistically sample from past reward ratios
        '''
        choice = np.random.rand() < self.prob
        self.choice_history.append(choice)
        return choice    
    

# A Matching agent object
class PiecewiseConstantProbAgent(Agent):
    '''
    Simulate an agent that decides with a fixed probability in blocks
    'Optimal' agent when it knows the probability of each block in a world
    '''
    def __init__(self, rates, ntrials):
        self.choice_history = []
        self.rates = rates
        self.ntrials = ntrials
        self.outcome_history = []
        self.curr_rate = rates[0]
        self.curr_block = 0
        
    def outcome_received(self, outcome):
        self.outcome_history.append(outcome)
        
    
    def make_choice(self):
        '''
        Make a choice, probabilistically sample from past reward ratios
        '''
        choice = np.random.rand() < self.curr_rate
        
        if len(self.outcome_history) > sum(self.ntrials[:self.curr_block + 1]):
            self.curr_block += 1
            self.curr_rate = self.rates[self.curr_block]
        
        self.choice_history.append(choice)
        return choice    

class EGreedyInferenceBasedAgent(Agent):
    def __init__(self, prew, pswitch, eps=0):
        '''
        An inference-based agent
        :param prew: probability of reward of high-state
        :param pswitch: probability of switch
        :param eps: value used for epsilon-greedy strategy
        '''
        self.prew = prew
        self.pswitch = pswitch
        self.choice_history = []
        self.p0_history = []
        self.p1_history = []
        self.eps = eps
        self.outcome_history = []
        self.Rewards0side_history = []
        self.Rewards1side_history = []
        self.type = type
        self.p0 = 0.5
        self.p1 = 0.5

    def outcome_received(self, outcome):
        self.p0, self.p1 = self.update_prob(self.choice_history[-1], outcome)

        self.outcome_history.append(outcome)
        self.p0_history.append(self.p0)
        self.p1_history.append(self.p1)

        # Update reward history on each side
        if outcome == 1:
            if self.choice_history[-1] == 1:
                self.Rewards1side_history.append(1)
                self.Rewards0side_history.append(0)
            else:
                self.Rewards0side_history.append(1)
                self.Rewards1side_history.append(0)


    def update_prob(self, choice, outcome):
        '''
        Returns the updated probability (p0, p1)
        :param choice: previous choice, 0 or 1
        :param outcome: previous outcome, 0 or 1
        :return: updated probability (p0, p1)
        '''
        prew = self.prew
        pswitch = self.pswitch

        if choice == 1 and outcome == 1:  # chose right, rew
            LLHrtGiven0 = 1 - prew
            LLHrtGiven1 = prew
        elif choice == 0 and outcome == 1:  # chose left, rew
            LLHrtGiven0 = prew
            LLHrtGiven1 = 1 - prew
        elif choice == 1 and outcome == 0:  # chose right, no rew
            LLHrtGiven0 = prew
            LLHrtGiven1 = 1 - prew
        else:  # chose left, no rew
            LLHrtGiven0 = 1 - prew
            LLHrtGiven1 = prew

        p1prev = self.p1
        p0prev = self.p0
        p0new = (1 - pswitch) * LLHrtGiven0 * p0prev + \
            pswitch * LLHrtGiven1 * p1prev
        p1new = pswitch * LLHrtGiven0 * p0prev + \
            (1 - pswitch) * LLHrtGiven1 * p1prev

        p0 = p0new / (p0new + p1new)
        p1 = p1new / (p0new + p1new)

        return p0, p1

    def make_choice(self):
        explore = np.random.rand() < self.eps
        if explore:
            choice = np.random.rand() < 0.5
        else:
            # Optimal agent picks the action with higher prob
            if self.p1 > self.p0:
                choice = 1
            elif self.p1 < self.p0:
                choice = 0
            else:
                choice = np.random.rand() > 0.5

            # Old code for random agent
            # choice = np.random.rand() > self.p0 / (self.p0 + self.p1)
        self.choice_history.append(choice)
        return choice
        
    
class Experiment():
    '''
    An experiment consists of an agent in a world
    '''
    def __init__(self, agent, world):
        self.agent = agent
        self.world = world
        
        
    def run(self):
        '''
        Run the proposed experiment
        '''
        choices = []
        rewards = []
        for i in range(sum(self.world.ntrials) - 1):
            choice = self.agent.make_choice()
            #print('choice = ', int(choice))
            choices.append(choice)
            reward = self.world.update(choice)
            #print('reward = ', int(reward))
            self.agent.outcome_received(reward)
            rewards.append(reward)
            
        return choices, rewards



    

class LocalMatchingAgent(Agent):
    '''
    Simulate an agent that matches with a leaky integrator
    '''
    def __init__(self, tau, eps):
        """
        tau: time constant of matching agent (integration kernel)
        eps: rate is limited to the range [eps, 1-eps]
        """
        assert(0 <= eps <= 0.5)
        self.tau = tau
        self.eps = eps
        self.Rewards1side = 1 # Start with 1 so that the agent does not get 'stuck'
        self.Rewards0side = 1
        self.choice_history = []
        self.Rewards1side_history = []
        self.Rewards0side_history = []
        self.outcome_history = []
        
    def outcome_received(self, outcome):
        if self.tau == 0:
            factor = 1
        else:
            factor = np.power(0.5, 1/self.tau)
        self.Rewards1side = self.Rewards1side * factor
        self.Rewards0side = self.Rewards0side * factor
        if outcome == 1:
            if self.choice_history[-1] == 1:
                self.Rewards1side += 1
            else:
                self.Rewards0side += 1
                                
        self.outcome_history.append(outcome)
        self.Rewards1side_history.append(self.Rewards1side)
        self.Rewards0side_history.append(self.Rewards0side)
              
    
    def make_choice(self):
        '''
        Make a choice, probabilistically sample from past reward ratios
        '''
        
        p = self.Rewards1side / (self.Rewards1side + self.Rewards0side)
        p = min(p, 1 - self.eps)
        p = max(p, self.eps)
        choice = np.random.rand() < p
        self.choice_history.append(choice)
        return choice
