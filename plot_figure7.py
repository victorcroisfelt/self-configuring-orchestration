import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Load data
data = np.load('data/figure7.npz')

n_pilot_subblocks = data['n_pilot_subblocks'] 
n_pilot_subblocks_probe_range = data['n_pilot_subblocks_probe_range'] / n_pilot_subblocks
eta_range = data['eta_range']
proba_false_alarm_range = data['proba_false_alarm_range']

pd_proba_detection = data['pd_proba_detection']
dsp_proba_detection = data['dsp_proba_detection']

#####

#
# Plot for eta
#

fig, ax = plt.subplots()

markers = ['s', 'd']
linestyles = ['-', '--']
colors = ['tab:orange', 'tab:green']

ax.plot(1 - eta_range, pd_proba_detection[7, :, 0]/100, linestyle=linestyles[0], color=colors[0], marker=markers[0], label=r'PD-enabled HRIS: $P_{\rm FA}=10^{-2}$')
ax.plot(1 - eta_range, dsp_proba_detection[7, :, 0]/100, linestyle=linestyles[0], color=colors[0], marker=markers[1], label=r'DSP-enabled HRIS: $P_{\rm FA}=10^{-2}$')

ax.plot(1 - eta_range, pd_proba_detection[7, :, 1]/100, linestyle=linestyles[1], color=colors[1], marker=markers[0], label=r'PD-enabled HRIS: $P_{\rm FA}=10^{-3}$')
ax.plot(1 - eta_range, dsp_proba_detection[7, :, 1]/100, linestyle=linestyles[1], color=colors[1], marker=markers[1], label=r'DSP-enabled HRIS: $P_{\rm FA}=10^{-3}$')

ax.set_xscale('log')

ax.set_xlabel(r'Fraction of power absorbed, $1-\eta$')
ax.set_ylabel(r'Probability of detection, $P_{\rm D}$')

ax.set_xscale('log')

#ax.set_ylim([0.55, 1.05])

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.legend(framealpha=0.5)

figure7a_pd_pfa_001 = np.vstack((1 - eta_range, pd_proba_detection[7, :, 0]/100)).T
np.savetxt('txts/fig7a_pd-hris_pfa001.txt', figure7a_pd_pfa_001, fmt='%.4f', comments='')

figure7a_pd_pfa_0001 = np.vstack((1 - eta_range, pd_proba_detection[7, :, 1]/100)).T
np.savetxt('txts/fig7a_pd-hris_pfa0001.txt', figure7a_pd_pfa_0001, fmt='%.4f', comments='')

figure7a_dsp_pfa_001 = np.vstack((1 - eta_range, dsp_proba_detection[7, :, 0]/100)).T
np.savetxt('txts/fig7a_dsp-hris_pfa001.txt', figure7a_dsp_pfa_001, fmt='%.4f', comments='')

figure7a_dsp_pfa_0001 = np.vstack((1 - eta_range, dsp_proba_detection[7, :, 1]/100)).T
np.savetxt('txts/fig7a_dsp-hris_pfa0001.txt', figure7a_dsp_pfa_0001, fmt='%.4f', comments='')

#
# Plot for relative probe duration
#

fig, ax = plt.subplots()

markers = ['s', 'd']
linestyles = ['-', '--']
colors = ['tab:orange', 'tab:green']

ax.plot(n_pilot_subblocks_probe_range, pd_proba_detection[:, -2, 0]/100, linestyle=linestyles[0], color=colors[0], marker=markers[0], label=r'PD-enabled HRIS: $P_{\rm FA}=10^{-2}$')
ax.plot(n_pilot_subblocks_probe_range, dsp_proba_detection[:, -2, 0]/100, linestyle=linestyles[0], color=colors[0], marker=markers[1], label=r'DSP-enabled HRIS: $P_{\rm FA}=10^{-2}$')

ax.plot(n_pilot_subblocks_probe_range, pd_proba_detection[:, -2, 1]/100, linestyle=linestyles[1], color=colors[1], marker=markers[0], label=r'PD-enabled HRIS: $P_{\rm FA}=10^{-3}$')
ax.plot(n_pilot_subblocks_probe_range, dsp_proba_detection[:, -2, 1]/100, linestyle=linestyles[1], color=colors[1], marker=markers[1], label=r'DSP-enabled HRIS: $P_{\rm FA}=10^{-3}$')

ax.set_xlabel(r'Relative probe duration, $\varrho$')
ax.set_ylabel(r'Probability of detection, $P_{\rm D}$')

#ax.set_ylim([0.55, 1.05])

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.legend(framealpha=0.5)

figure7b_pd_pfa_001 = np.vstack((n_pilot_subblocks_probe_range, pd_proba_detection[:, -2, 0]/100)).T
np.savetxt('txts/fig7b_pd-hris_pfa001.txt', figure7b_pd_pfa_001, fmt='%.4f', comments='')

figure7b_pd_pfa_0001 = np.vstack((n_pilot_subblocks_probe_range, pd_proba_detection[:, -2, 1]/100)).T
np.savetxt('txts/fig7b_pd-hris_pfa0001.txt', figure7b_pd_pfa_0001, fmt='%.4f', comments='')

figure7b_dsp_pfa_001 = np.vstack((n_pilot_subblocks_probe_range, dsp_proba_detection[:, -2, 0]/100)).T
np.savetxt('txts/fig7b_dsp-hris_pfa001.txt', figure7b_dsp_pfa_001, fmt='%.4f', comments='')

figure7b_dsp_pfa_0001 = np.vstack((n_pilot_subblocks_probe_range, dsp_proba_detection[:, -2, 1]/100)).T
np.savetxt('txts/fig7b_dsp-hris_pfa0001.txt', figure7b_dsp_pfa_0001, fmt='%.4f', comments='')

print('--------------------')
print('Probability of Detection')
print('--------------------')

print('PD-enabled = {:.2f}'.format(pd_proba_detection[7, -2, 0]))
print('DSP-enabled = {:.2f}'.format(dsp_proba_detection[7, -2, 0]))

plt.show()