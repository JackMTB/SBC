# Model of SBC

## Theory

The following investigation aims at determing the pulse time and number of pulses which leads to the highest cooling rate in SBC.
Sideband cooling works by driving in cycles the red sideband transition, which excites the transition

$ |g,n\rangle \leftrightarrow |e,n-1\rangle $

such that $\bar{n}$ reduces at each cycle. After Doppler cooling, the ion is typically in a mixed state, which can be described by

$$ | \Psi_{DC}\rangle = \sum_{n} p_{\bar{n}}(n) |g,n\rangle \langle g,n|, $$

where

$$p_{\bar{n}}(n) = \frac{\bar{n}^n}{(\bar{n} + 1)^{n+1}}$$

is the probability of being in a Fock state $n$ given a thermal distribution characterised by $\bar{n}$.

The Rabi frequency of red sideband transition is dependent on the fock state

$$ \Omega_{n,n-1} = \langle n | e^{i \eta (\hat{a} + \hat{a}^\dagger)} |n-1\rangle. $$

This implies that if we drive RSB transistion each state will oscillate at its respective rabi frequency, therefore the signal observed is given by a superposition of oscillations at different frequencies, each weighted by the respective probability densisty. Assuming a thermal distribution, the signal measured upon driving the RSB on resonance is going to be

$$ p_{|e>}(t) = \sum_{n} p_{\bar{n}}(n)  p_{exc}(t,n),$$

where

$$ p_{exc}(t,n) = \sin^2\left( \frac{\Omega_{n,n-1} t}{2} \right)$$

is the probability of excitation from state $|g,n\rangle \leftrightarrow |e,n-1\rangle$ for a given pulse time $t$.

## Model for simulating evolution of initial distribution given pulse sequence

The model takes as input a pulse sequence, an initial $\bar{n}$ and $\Omega_0$ and will output a matrix where each row is the distribution of fock states after applying a pulse of the sequence.

for a pulse $i$ applied for a time $t$, the evolution of the population in the fock state $n$ is

$$ pop(i,n) = pop(i-1,n) + pop(i-1,n+1) \cdot p_{exc}(t,n+1) - pop(i-1,n) \cdot p_{exc}(t,n).$$

by specifying the argument simulation_type = 'monte carlo', the same concept is implemented with a monte carlo approach. However, this method is more computationally expensive.