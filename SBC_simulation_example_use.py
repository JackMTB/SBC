import SBC_simulation as sbc
import numpy as np
import matplotlib.pyplot as plt

# Initial parameters
rabi_fr = 2 * np.pi * 100e3
initial_nbar = 30
LD_param = 0.1

# Create a sideband cooling sequence
seq_jj = sbc.create_sequence_n1(initial_nbar, rabi_fr, fidelity=1e-3)

# initialise figure
pulse_sequence_list = [seq_jj,seq_jj]
labels = ['motional decay dispersion: OFF','motional decay dispersion: ON']
col = ['#EF7C8E','#2E8BC0']

fig = plt.figure(figsize=(15,25))
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
#-----------------------------------------------------------
# perform two different types of simulations:
i=0
for pulse_sequence in pulse_sequence_list:
    #-------------------------------------------------------
    # simulation 
    sequence_simulation = sbc.Simulate_sequence(
                pulse_sequence,
                initial_nbar,
                rabi_fr,
                LD_param,
                upper_sb_transition = 2,
                dispersion=i)
    #-------------------------------------------------------
    # Analysis
    analysis = sbc.Analyse_sequence(sequence_simulation.distribution)
    #-------------------------------------------------------
    # Plot
    ax1.plot(analysis.nbar,linewidth=2, color=col[i], label = labels[i])
    ax2.bar(np.arange(1,len(analysis.final_distribution),1)+i*0.4, analysis.final_distribution[1:], width = 0.3, color=col[i])
    ax3.plot(1-analysis.p_0,linewidth=2, color=col[i])
    ax4.plot(analysis.variance,linewidth=2, color=col[i])
    i+=1


ax1.set_xlabel('Pulse number', fontsize=18)
ax1.set_ylabel(r'Temperature $\bar{n}$', fontsize=18)
ax1.tick_params(axis='both', labelsize= 18)
ax1.set_xlim(0)
ax1.set_yscale('log')
ax1.legend(fontsize=18)
ax1.grid()

ax2.set_xlabel(r'$n$', fontsize=18)
ax2.set_ylabel('P(n)', fontsize=18)
ax2.tick_params(axis='both', labelsize= 18)
ax2.set_ylim(0)
ax2.set_title('Population starting from n=1 at end of sequence', fontsize = 18)
ax2.grid()

ax3.set_xlabel('Pulse number', fontsize=18)
ax3.set_ylabel('1-Population in n=0', fontsize=18)
ax3.tick_params(axis='both', labelsize= 18)
ax3.set_xlim(0)
ax3.set_yscale('log')
ax3.grid()
ax4.set_xlabel('Pulse number', fontsize=18)
ax4.set_ylabel('Variance wrt\n thermal distribution', fontsize=18)
ax4.tick_params(axis='both', labelsize= 18)
ax4.set_xlim(0)
ax4.set_ylim(0)
ax4.grid()

plt.show()