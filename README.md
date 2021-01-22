# Model of SBC

## Theory

The following investigation aims at determing the pulse time and number of pulses which leads to the highest cooling rate in SBC.
Sideband cooling works by driving in cycles the red sideband transition, which excites the transition

$$ |g,n\rangle \leftrightarrow |e,n-1\rangle $$

such that $\bar{n}$ reduces at each cycle. After Doppler cooling, the ion is typically in a mixed state, which can be described by

![eq2](https://latex.codecogs.com/gif.latex?%7C%20%5CPsi_%7BDC%7D%5Crangle%20%3D%20%5Csum_%7Bn%7D%20p_%7B%5Cbar%7Bn%7D%7D%28n%29%20%7Cg%2Cn%5Crangle%20%5Clangle%20g%2Cn%7C%2C)

where

![eq3](https://latex.codecogs.com/gif.latex?p_%7B%5Cbar%7Bn%7D%7D%28n%29%20%3D%20%5Cfrac%7B%5Cbar%7Bn%7D%5En%7D%7B%28%5Cbar%7Bn%7D%20&plus;%201%29%5E%7Bn&plus;1%7D%7D)

is the probability of being in a Fock state $n$ given a thermal distribution characterised by $\bar{n}$.

The Rabi frequency of red sideband transition is dependent on the fock state

![eq4](https://latex.codecogs.com/gif.latex?%5COmega_%7Bn%2Cn-1%7D%20%3D%20%5Clangle%20n%20%7C%20e%5E%7Bi%20%5Ceta%20%28%5Chat%7Ba%7D%20&plus;%20%5Chat%7Ba%7D%5E%5Cdagger%29%7D%20%7Cn-1%5Crangle.)

This implies that if we drive RSB transistion each state will oscillate at its respective rabi frequency, therefore the signal observed is given by a superposition of oscillations at different frequencies, each weighted by the respective probability densisty. Assuming a thermal distribution, the signal measured upon driving the RSB on resonance is going to be

![eq5](https://latex.codecogs.com/gif.latex?p_%7B%7Ce%3E%7D%28t%29%20%3D%20%5Csum_%7Bn%7D%20p_%7B%5Cbar%7Bn%7D%7D%28n%29%20p_%7Bexc%7D%28t%2Cn%29%2C)

where

![eq6](https://latex.codecogs.com/gif.latex?p_%7Bexc%7D%28t%2Cn%29%20%3D%20%5Csin%5E2%5Cleft%28%20%5Cfrac%7B%5COmega_%7Bn%2Cn-1%7D%20t%7D%7B2%7D%20%5Cright%29)

is the probability of excitation from state $|g,n\rangle \leftrightarrow |e,n-1\rangle$ for a given pulse time $t$.

## Model for simulating evolution of initial distribution given pulse sequence

The model takes as input a pulse sequence, an initial $\bar{n}$ and $\Omega_0$ and will output a matrix where each row is the distribution of fock states after applying a pulse of the sequence.

for a pulse $i$ applied for a time $t$, the evolution of the population in the fock state $n$ is

![eq7](https://latex.codecogs.com/gif.latex?pop%28i%2Cn%29%20%3D%20pop%28i-1%2Cn%29%20&plus;%20pop%28i-1%2Cn&plus;1%29%20%5Ccdot%20p_%7Bexc%7D%28t%2Cn&plus;1%29%20-%20pop%28i-1%2Cn%29%20%5Ccdot%20p_%7Bexc%7D%28t%2Cn%29.)

by specifying the argument simulation_type = 'monte carlo', the same concept is implemented with a monte carlo approach. However, this method is more computationally expensive.
