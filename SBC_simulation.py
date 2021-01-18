import numpy as np
import scipy.optimize as opt
import scipy.special as sp

def product_pos_m(n, m):
    if type(n)==int:
        prod = 1
    else:
        prod = np.ones(np.shape(n), dtype=int)
    for i  in range(0,m):
        prod = prod / (n+m-i)
    return prod**0.5

def product_neg_m(n, m):
    if type(n)==int:
        prod = 1
        if n-m>=0:
            for i  in range(0,m):
                prod = prod / (n+m-i)
        else:
            prod = 0
    else:
        prod = np.ones(np.shape(n), dtype=float)
        for i  in range(0,m):
            prod[n-m>=0] = prod[n-m>=0] / (n[n-m>=0]+m-i)
        prod[n-m<0] = 0
        

    # else:
    #     prod = np.ones(np.shape(n), dtype=float)
    #     cacca = np.ones(np.shape(n[n-m>=0]),dtype=float)
    #     for i  in range(0,m):
    #         cacca = cacca / (n[n-m>=0]+2*m-i)
    #     prod[n-m>=0] = prod[n-m>=0] * cacca
    #     prod[n-m<0] = 0

    return prod**0.5

def Om_n_m(eta,n,m,Om0):

    if type(n)==int:
        if n<0:
            output = 0

        else:
            t1 = np.exp(-eta**2/2)
            t2 = eta**np.abs(m)
            if m == 0:
                t3 = 1
                t4 = np.abs(sp.assoc_laguerre(eta**2,n,np.abs(m)))
                
            elif m > 0:
                t3 = product_pos_m(n, m)
                t4 = np.abs(sp.assoc_laguerre(eta**2,n,np.abs(m)))
            
            elif (m < 0) and (n+m>=0):
                t3 = product_neg_m(n, -m)
                t4 = np.abs(sp.assoc_laguerre(eta**2,n,np.abs(m)))

            elif (m < 0) and (n+m<0):
                t3 = 1
                t4 = 0

            output = Om0*t1*t2*t3*t4
    else:
        output = np.zeros(len(n))

        t1 = np.exp(-eta**2/2)
        t2 = eta**np.abs(m)
        if m == 0:
            t3 = 1
            t4 = np.abs(sp.assoc_laguerre(eta**2,n[n>=0],np.abs(m)))
            
        elif m > 0:
            t3 = product_pos_m(n[n>=0], m)
            t4 = np.abs(sp.assoc_laguerre(eta**2,n[n>=0],np.abs(m)))
        
        elif m < 0:
            t3 = product_neg_m(n[n>=0], -m)
            t4 = np.abs(sp.assoc_laguerre(eta**2,n[n>=0],np.abs(m)))

        output[n>=0] = Om0*t1*t2*t3*t4

    return output

class prob_motion_decay:

    def __init__(self,n, LD_param, upper_sb_transition = 7):
        self.upper_sb_transition = upper_sb_transition
        self.n = n
        #list of m-th addressed sidebands
        if upper_sb_transition==0:
            m_list = [0]
        else:
            m_list = np.arange(-upper_sb_transition,upper_sb_transition+1, dtype = int)

        #initialise matrix
        all_probs = np.zeros((len(m_list),len(n)))

        #calculate the modified_Rabi given a list of fock state "n", for a given m, for all m.
        for i, m in enumerate(m_list):
            all_probs[i] = sbc.Om_n_m(LD_param, n+m, -m, 1)
        
        #normalise the probability. sum of p over all m =1
        normalisation = sum(all_probs)
        all_probs[:,n>0] = all_probs[:,n>0] / normalisation[n>0]

        #assign the probability matrix
        self.all_probs = all_probs

    def delta_n(self,sb):
        #return the probability for the list of fock state given
        #  a transition to the m-th sb
        if (np.shape(self.all_probs)[0] == 1) and (sb!=0):
            output = np.zeros(len(self.n))
        else:
            output = self.all_probs[self.upper_sb_transition + sb]
        return output

def therm_dist(n,nb):
    """
    Return occupation probability of state n for thermal dist, nbar.
    """
    therm_ratio = nb/(nb+1)
    a_n = therm_ratio**n/(nb+1)
    return a_n

def get_a_param(rabi_fr):
    return -1.7670938755403408e-13 + 7.828211198367295 * rabi_fr**(-0.9999999948137507)

def get_b_param(rabi_fr):
    return 2.1187600574881472e-13 + 21.160046203613323 * rabi_fr**(-1.0000000020859856)

def get_optimal_pi_time(nbar,rabi_fr):
    a = get_a_param(rabi_fr)
    b = get_b_param(rabi_fr)
    c = 0.38674354614988055
    
    return a + b*np.exp(-c*nbar)

def calculate_upper_lim1(nbar, fidelity):
        n=0
        prob_density_of_n = 2
        while prob_density_of_n > fidelity:
            prob_density_of_n = therm_dist(n,nbar)
            n+=1
        return n

def calculate_upper_n(nbar, fidelity):
        n = 1
        reached_fidelity = 1
        while reached_fidelity > fidelity:
            n_list = np.arange(0,n,1)
            prob_density_of_n = therm_dist(n_list,nbar)
            reached_fidelity = 1-sum(prob_density_of_n)
            n += 1
        return n

def generate_sequence(weights,rabi_fr):
    sequence = np.zeros(int(sum(weights)))
    i = 0
    for n, weight in enumerate(weights):

        if weight!=0:
            sequence[i:i+weight] = np.full(weight, get_optimal_t(n+1, rabi_fr))
            i+= weight
    return sequence[::-1]

def create_sequence_n1(nbar, rabi_fr, fidelity = 1e-4):
    weights = np.ones(calculate_upper_n(nbar, fidelity), dtype=int)
    print(len(weights))
    sequence = generate_sequence(weights,rabi_fr)
    return sequence

def get_optimal_t(n, rabi_fr, LD_param=0.107):
    return np.pi / Om_n_m(LD_param, n, -1, rabi_fr)

def prob_excitation_n(t, n_list, rabi_fr, LD_param):
        rabi_rsb = Om_n_m(LD_param, n_list, -1, rabi_fr)
        return np.sin(rabi_rsb*0.5*t)**2

def create_const_seq(nbar, rabi_fr, n_pulses, fraction = 0.95):

    seq = np.zeros(n_pulses)
    seq[:int(fraction * n_pulses)] = np.full(int(fraction * n_pulses), get_optimal_pi_time(nbar,rabi_fr))
    seq[int(fraction * n_pulses):] = np.full(n_pulses-int(fraction * n_pulses), get_optimal_t(1, rabi_fr))
    return seq

def create_lin_seq(nbar, rabi_fr, n_pulses):
    min_t = get_optimal_pi_time(nbar,rabi_fr)
    max_t = get_optimal_t(1, rabi_fr)
    seq = np.linspace(min_t, max_t, n_pulses)
    return seq



class Simulate_sequence:
    # Simulates how a pulse sequence impacts the initial thermal distribution of an ion
    
    def __init__(self,pulse_sequence, nbar, rabi_fr, LD_param,
                upper_sb_transition = 6,
                simulation_type='fast',
                dispersion=False,
                initial_distribution=None,
                max_fidelity=1e-6): #------------------------------------------------------------change back
        self.dispersion = dispersion
        self.upper_sb_transition = upper_sb_transition
        self.LD_param = LD_param
        self.nbar = nbar
        self.rabi_fr = rabi_fr
        self.pulse_sequence = pulse_sequence
        self.n_pulses_4_pn01 = 0
        self.max_fidelity = max_fidelity # corresponds to 1-sum(thermal pop considered)
        
        self.initial_distribution = initial_distribution
        if simulation_type == 'fast':
            self.simulate_cooling_sequence()
        elif simulation_type == 'monte carlo':
            self.simulate_cooling_sequence_MC()
        
    def simulate_transition_occurence(self, prob_excitation):
        # returns an array where each element represents whether a transition from n to n-1 has occurred
        simulated_outcome = np.random.rand(len(prob_excitation))
        transitions_occured = np.zeros(len(prob_excitation),dtype=bool)
        transitions_occured[simulated_outcome<=prob_excitation] = True
        return transitions_occured

    def update_population(self, current_thermal_dist, occured_transitions):
        for i in range(len(occured_transitions)):
            if occured_transitions[i]==True:
                current_thermal_dist[i-1] += current_thermal_dist[i]
                current_thermal_dist[i] = 0
        return current_thermal_dist

    def simulate_cooling_sequence_MC(self):
        
        n = np.arange(0,calculate_upper_n(self.nbar, self.max_fidelity),1)
        distribution = np.zeros((len(self.pulse_sequence)+1,len(n)))
        summed_distribution = np.zeros((len(self.pulse_sequence)+1,len(n)))
        
        num_repeat = 100
        for i in range(num_repeat):
            
            if self.initial_distribution == None:
                distribution[0] = therm_dist(n,nbar)
            else:
                distribution[0] = self.initial_distribution
                
            for i, pulse_time in enumerate(self.pulse_sequence):
                prob_exc = prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)
                occurred_transitions = self.simulate_transition_occurence(prob_exc)
                distribution[i+1] = self.update_population(distribution[i], occurred_transitions)
            
            summed_distribution += distribution
            
        self.distribution = summed_distribution / num_repeat
    

        
    def simulate_cooling_sequence(self):
        
        # create a list representing the fock states taken into consideration
        n = np.arange(0,calculate_upper_n(self.nbar, self.max_fidelity),1)

        # create a list that contains the total prob of excitation that each fock
        # state received throughout the pulse sequence. For analysis only, not
        # actually used for the simulation.
        self.accumulated_prob = np.zeros(len(n))

        # create a matrix that will contain the evolution of the population distribution throughout the pulse sequence
        distribution = np.zeros((len(self.pulse_sequence)+1,len(n)))

        # Initialise the initial distribution, thermal by default.
        if self.initial_distribution == None:
            distribution[0] = therm_dist(n,self.nbar)
        else:
            distribution[0] = self.initial_distribution

        # Calculate the prob of decay to the mth sideband
        prob_motion_D = prob_motion_decay(n, self.LD_param, upper_sb_transition= self.upper_sb_transition)
        
        # Calculate variation in population distribution for each pulse in the pulse sequence.
        for i, pulse_time in enumerate(self.pulse_sequence):

            prob_exc = prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)
            self.accumulated_prob += prob_exc 
            
            if self.dispersion == False:
                positive_contribution = np.zeros(len(n))
                positive_contribution[:-1] = distribution[i,1:] * prob_excitation_n(pulse_time, n[1:], self.rabi_fr, self.LD_param)
                negative_contribution = distribution[i] * prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)

                distribution[i+1] = distribution[i] + positive_contribution - negative_contribution

            # elif self.dispersion == True:
            #     extended_pop = np.zeros(len(n)+2*self.upper_sb_transition+1)
            #     extended_pop[self.upper_sb_transition+1:-self.upper_sb_transition] = distribution[i]
                
            #     extended_exc_prob = np.zeros(len(n)+2*self.upper_sb_transition+1)
            #     extended_exc_prob[self.upper_sb_transition+1:-self.upper_sb_transition] = prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)

            #     positive_contribution = np.zeros(len(n))
            #     extended_mot_prob = np.zeros(len(n)+2*self.upper_sb_transition)
            #     for m in np.arange(-self.upper_sb_transition,self.upper_sb_transition+1,dtype=int):
                    
            #         extended_mot_prob[self.upper_sb_transition+1:-self.upper_sb_transition] = prob_motion_D.delta_n(m)
            #         print(f'')
            #         positive_contribution = positive_contribution +\
            #                                 extended_pop[(self.upper_sb_transition-m+1) : (self.upper_sb_transition+len(n)-m+1)] *\
            #                                 extended_exc_prob[(self.upper_sb_transition-m+1) : (self.upper_sb_transition+len(n)-m+1)] *\
            #                                 extended_mot_prob[(self.upper_sb_transition-m) : (self.upper_sb_transition+len(n)-m)]

            #     negative_contribution = np.zeros(len(n))
            #     negative_contribution[1:] = distribution[i,1:] * prob_excitation_n(pulse_time, n[1:], self.rabi_fr, self.LD_param) * (1-prob_motion_D.delta_n(1)[:-1])
                
            #     print(sum(positive_contribution)/sum(negative_contribution))

            #     distribution[i+1] = distribution[i] + positive_contribution - negative_contribution

            elif self.dispersion == True:
                pos_contribution = np.zeros(len(n))
                neg_contribution = np.zeros(len(n))
                prob_exc = prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)

                for ni in n:
                    
                    for m in np.arange(-self.upper_sb_transition,self.upper_sb_transition+1,dtype=int):
                        prob_mot = prob_motion_D.delta_n(m)

                        if (ni+m+1<0) or (ni+m+1>len(n)-1):
                            continue
                        else:
                            pos_contribution[ni] +=  prob_exc[ni+m+1] * distribution[i,ni+m+1] * prob_mot[ni]

                    
                    prob_mot = prob_motion_D.delta_n(1)
                    neg_contribution[ni] = prob_exc[ni] * distribution[i,ni] * (1-prob_mot[ni])

                print(sum(pos_contribution)/sum(neg_contribution))              
                distribution[i+1] = distribution[i] + pos_contribution - neg_contribution

        self.distribution = distribution

class Analyse_sequence:
    
    def __init__(self,distribution):
        self.number_n_included = len(distribution[0])
        self.number_pulses = len(distribution)
        self.pop_included = sum(distribution[0])
        
        self.analyse_distribution(distribution)
        self.final_distribution = distribution[-1]
        
        
    def get_nbar(self, t, distribution):
        return sum(distribution[t] * np.arange(0,self.number_n_included,1))
        
    def get_p0(self, t, distribution):
        return distribution[t,0]
    
    def variance_wrt_thermal(self, t, nbar, distribution):
        data = distribution[t]
        thermal_model = therm_dist(np.arange(0,self.number_n_included,1), nbar)
        
        variance = np.zeros(len(data))
        for i in range(len(data)):
            variance[i] = (data[i]-thermal_model[i])**2
        return sum(variance)/len(variance)

    def analyse_distribution(self,distribution):

        nbar = np.zeros(self.number_pulses)
        p_0 = np.zeros(self.number_pulses)
        variance = np.zeros(self.number_pulses)

        for t in range(self.number_pulses):
            nbar[t] = self.get_nbar(t, distribution)
            p_0[t] = self.get_p0(t, distribution)
            variance[t] = self.variance_wrt_thermal(t, nbar[t], distribution)

        self.nbar = nbar
        self.p_0 = p_0
        self.variance = variance