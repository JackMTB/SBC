import numpy as np
import scipy.optimize as opt
import scipy.special as sp

#----------------Modified Rabi------------------------------
def Om_n_m(eta,n,m,Om0):
    # Calculates the modified Rabi frequency associated with the transition between |g,n> -- |e,n+m>
    # accepts both integer and np.arrays
    # also work for n+m<0
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
            output[n>=0] = Om0*t1*t2*t3*t4
            
        elif m > 0:
            t3 = product_pos_m(n[n>=0], m)
            t4 = np.abs(sp.assoc_laguerre(eta**2,n[n>=0],np.abs(m)))
            output[n>=0] = Om0*t1*t2*t3*t4
        
        elif m < 0:
            t3 = product_neg_m(n[n>=0], -m)
            t4 = np.abs(sp.assoc_laguerre(eta**2,n[n>=0],np.abs(m)))
            output_temp = np.zeros(len(n))
            output_temp[n>=0] = Om0*t1*t2*t3*t4
            output[n+m>=0] = output_temp[:len(output[n+m>=0])]
    return output

def product_pos_m(n, m):
    # used to calculate part of the modified rabi frequency for positive m
    if type(n)==int:
        prod = 1
    else:
        prod = np.ones(np.shape(n), dtype=int)
    for i  in range(0,m):
        prod = prod / (n+m-i)
    return prod**0.5

def product_neg_m(n, m):
    # used to calculate part of the modified rabi frequency for negative m
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
            prod = prod / (n+m-i)
    return prod**0.5

    
#-------------Probability of motional decay-----------------
class prob_motion_decay:
    # Create an object containing P_mot(n,m) the probability of DECAY from |e,n> -- |g,n+m>.
    # Given a maximum number of considered sideband transitions -M-, 
    # the probability is normalised as sum_{m=-M}^{m=M} P_mot(n,m) = 1.
    # call as:
    #   probability_motion = prob_motion_decay(n, LD_param)
    #   probability_motion_m = probability_motion.delta_n(m)
    # accepts n as integer or numpy array
    def __init__(self,n, LD_param, upper_sb_transition = 7):
        self.upper_sb_transition = upper_sb_transition
        self.n = n
        self.LD_param = LD_param
        #list of m-th addressed sidebands
        if upper_sb_transition==0:
            m_list = [0]
        else:
            m_list = np.arange(-upper_sb_transition,upper_sb_transition+1, dtype = int)

        self.normalisation = np.zeros(len(n))
        for m in m_list:
            self.normalisation = self.normalisation + Om_n_m(LD_param,n-m,m,1)

    def delta_n(self,sb):
        if (self.upper_sb_transition == 0) and (sb!=0):
            output = np.zeros(len(self.n))
        else:
            output = np.zeros(len(self.n))
            output[self.n>=0] = Om_n_m(self.LD_param,self.n[self.n>=0]+sb,-sb,1) / self.normalisation[self.n>=0]
        return output


#-----------------------------------------------------------
# Function for generation of SBC sequences and auxilliary 
# functions of the simulation
#-----------------------------------------------------------

def therm_dist(n,nb):
    # Return occupation probability of state n for a thermal distribution characterised by nbar.
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

def calculate_upper_n(nbar, fidelity):
    # calculates the N for which sum_{n=0,N}(thermal_dist(n)) = 1-fidelity
        n = 1
        reached_fidelity = 1
        while reached_fidelity > fidelity:
            n_list = np.arange(0,n,1)
            prob_density_of_n = therm_dist(n_list,nbar)
            reached_fidelity = 1-sum(prob_density_of_n)
            n += 1
        return n

def generate_sequence(weights,rabi_fr,LD_param):
    sequence = np.zeros(int(sum(weights)))
    i = 0
    for n, weight in enumerate(weights):

        if weight!=0:
            sequence[i:i+weight] = np.full(weight, get_optimal_t(n+1, rabi_fr,LD_param))
            i+= weight
    return sequence[::-1]

def create_sequence_n1(LD_param, rabi_fr, upper_n, sb=-1):
    list_n = np.arange(1,upper_n,dtype=int)
    sequence = np.pi / Om_n_m(LD_param, list_n,sb, rabi_fr)
    sb_order = np.ones(len(list_n)) * (-1)
    return [sequence[::-1],sb_order]

def get_optimal_t(n, rabi_fr, LD_param):
    return np.pi / Om_n_m(LD_param, n, -1, rabi_fr)

def prob_excitation_n(t, n_list, rabi_fr, LD_param, sb=-1):
        rabi_rsb = Om_n_m(LD_param, n_list, sb, rabi_fr)
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

#---------------Actual simulation---------------------------
class Simulate_sequence:
    # Simulates how a pulse sequence impacts the initial thermal distribution of an ion
    
    def __init__(self,pulse_sequence, nbar, rabi_fr, LD_param,
                upper_sb_transition = 6,
                simulation_type='fast',
                dispersion=False,
                initial_distribution=None,
                max_fidelity=1e-6): 
        self.dispersion = dispersion
        self.upper_sb_transition = upper_sb_transition
        self.LD_param = LD_param
        self.nbar = nbar
        self.rabi_fr = rabi_fr
        self.pulse_sequence = pulse_sequence[0]
        self.sb_order = pulse_sequence[1] 
        self.n_pulses_4_pn01 = 0
        self.max_fidelity = max_fidelity # corresponds to 1-sum(thermal pop considered)
        
        self.initial_distribution = initial_distribution
        if simulation_type == 'fast':
            self.simulate_cooling_sequence()
        elif simulation_type == 'monte carlo':
            self.simulate_cooling_sequence_MC()
        
    def simulate_transition_occurence(self, prob_excitation):
        # ONLY FOR MONTE CARLO
        # returns an array where each element represents whether a transition from n to n-1 has occurred
        simulated_outcome = np.random.rand(len(prob_excitation))
        transitions_occured = np.zeros(len(prob_excitation),dtype=bool)
        transitions_occured[simulated_outcome<=prob_excitation] = True
        return transitions_occured

    def update_population(self, current_thermal_dist, occured_transitions):
        # ONLY FOR MONTE CARLO
        for i in range(len(occured_transitions)):
            if occured_transitions[i]==True:
                current_thermal_dist[i-1] += current_thermal_dist[i]
                current_thermal_dist[i] = 0
        return current_thermal_dist

    def simulate_cooling_sequence_MC(self):
        # Monte carlo type simulation
        # not updated, not working - do not use
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
        extra_pop=0
        for i, pulse_time in enumerate(self.pulse_sequence):

            chosen_exc_sb = self.sb_order[i]
            prob_exc = prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)
            self.accumulated_prob += prob_exc 
            
            if self.dispersion == False:
                
                if chosen_exc_sb == -1:

                    positive_contribution = np.zeros(len(n))
                    positive_contribution[:-1] = distribution[i,1:] * prob_excitation_n(pulse_time, n[1:], self.rabi_fr, self.LD_param)
                    negative_contribution = distribution[i] * prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)
                
                elif chosen_exc_sb == -2:

                    positive_contribution = np.zeros(len(n))
                    positive_contribution[:-2] = distribution[i,2:] * prob_excitation_n(pulse_time, n[2:], self.rabi_fr, self.LD_param, -2)
                    negative_contribution = distribution[i] * prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param, -2)
               
                distribution[i+1] = distribution[i] + positive_contribution - negative_contribution

            elif self.dispersion == True:

                if self.upper_sb_transition>=1:
                    extended_pop = np.zeros(len(n)+2*self.upper_sb_transition)
                    extended_pop[self.upper_sb_transition-1 : self.upper_sb_transition-1+len(n)] = distribution[i]
                
                    extended_exc_prob = np.zeros(len(n)+2*self.upper_sb_transition)
                    extended_exc_prob[self.upper_sb_transition-1 : self.upper_sb_transition-1+len(n)] = prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param)

                    m_list = np.arange(-self.upper_sb_transition,self.upper_sb_transition+1,dtype=int)
                    m_list = np.delete(m_list,self.upper_sb_transition+1)

                elif self.upper_sb_transition==0:
                    extended_pop = np.zeros(len(n))
                    extended_pop[:-1] = distribution[i,1:]
                
                    extended_exc_prob = np.zeros(len(n))
                    extended_exc_prob[:-1] = prob_excitation_n(pulse_time, n[1:], self.rabi_fr, self.LD_param)

                    m_list = [0]

                positive_contribution = np.zeros(len(n))
                extended_mot_prob = np.zeros(len(n)+2*self.upper_sb_transition)
                for m in m_list:

                    extended_mot_prob[self.upper_sb_transition : self.upper_sb_transition + len(n)] = prob_motion_D.delta_n(m)

                    positive_contribution = positive_contribution +\
                                            extended_pop[self.upper_sb_transition-m : self.upper_sb_transition+len(n)-m] *\
                                            extended_exc_prob[self.upper_sb_transition-m : self.upper_sb_transition+len(n)-m] *\
                                            extended_mot_prob[self.upper_sb_transition-m : self.upper_sb_transition+len(n)-m]

                negative_contribution = np.zeros(len(n))
                extended_mot_prob = np.zeros(len(n)+1)
                extended_mot_prob[1:] = prob_motion_D.delta_n(1)
                negative_contribution = distribution[i] * prob_excitation_n(pulse_time, n, self.rabi_fr, self.LD_param) * (1-extended_mot_prob[:-1])

                distribution[i+1] = distribution[i] + positive_contribution - negative_contribution

        self.distribution = distribution

#------------Analysis of the simulation---------------------
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
