import SBC_simulation as sbc
import numpy as np
import matplotlib.pyplot as plt


rabi_fr = 2 * np.pi * 100e3
initial_nbar = 30

#my sequence
seq_jj = sbc.create_sequence_n1(initial_nbar, rabi_fr, fidelity=1e-3)
seq_const = sbc.create_const_seq(initial_nbar, rabi_fr, len(seq_jj))
seq_const1 = sbc.create_lin_seq(initial_nbar+3, rabi_fr, len(seq_jj))
seq_const2 = sbc.create_lin_seq(initial_nbar-5, rabi_fr, len(seq_jj))

pulse_sequence_list = [seq_jj, seq_const, seq_const1, seq_const2]
# nbars=[initial_nbar,initial_nbar-10,initial_nbar+10]
labels = ['my seq', 'constant seq', 'c+3','c-3']

fig = plt.figure(figsize=(15,7))
ax1 = fig.add_subplot(111)
# ax2 = fig.add_subplot(212)
i=0
for pulse_sequence in pulse_sequence_list:
    # simulation and analysis
    sequence_simulation = sbc.Simulate_sequence(pulse_sequence, initial_nbar, rabi_fr)
    analysis = sbc.Analyse_sequence(sequence_simulation.distribution)

    # plot
    ax1.plot(1-analysis.p_0,linewidth=2, label = labels[i])
    # ax1.plot(pulse_sequence*1e6,linewidth=2, label = labels[i])
    # ax2.bar(np.arange(1,len(analysis.final_distribution),1)+i*0.5, analysis.final_distribution[1:], width = 0.3)
    # note: the initial nbar doesn't correspond to the input nbar because of the truncation of the thermal distribution
    print('Fidelity of ground state preparation:' + '{:.3e}'.format(1-analysis.p_0[-1]))
    print('Final nbar:' + '{:.3e}'.format(analysis.nbar[-1]))
    print('Number of pulses in sequence: ' + str(len(pulse_sequence)))
    i+=1

ax1.set_xlabel('Pulse number', fontsize=18)
ax1.set_ylabel('1-Population in n=0', fontsize=18)

# ax1.set_xlabel('Pulse number', fontsize=18)
# ax1.set_ylabel('Time [us]', fontsize=18)

ax1.tick_params(axis='both', labelsize= 18)
# ax1.set_ylim(top=1)
ax1.set_xlim(0)
ax1.set_yscale('log')
ax1.grid()
ax1.legend(fontsize = 18)

plt.show()