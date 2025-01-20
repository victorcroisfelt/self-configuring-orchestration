import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Load data
data = np.load('data/figure7.npz')

n_pilot_subblocks_probe = data['n_pilot_subblocks_probe'] 
n_pilot_subblocks_probe_range = data['n_pilot_subblocks_probe_range'] / n_pilot_subblocks_probe
eta_range = data['eta_range']
proba_false_alarm_range = data['proba_false_alarm_range']

pow_pd = data['pow_pd']
sig_pd = data['sig_pd']

pow_nmse = data['pow_nmse2']
sig_nmse = data['sig_nmse2']

# Get average values
avg_pow_pd = pow_pd.mean(axis=(-1, -2))
avg_sig_pd = sig_pd.mean(axis=(-1, -2))

avg_pow_nmse = pow_nmse.mean(axis=(-1, -2))
avg_sig_nmse = sig_nmse.mean(axis=(-1, -2))

# Get worst case probability of detection
worst_pow_pd = np.nanmin(avg_pow_pd, axis=2)
worst_sig_pd = np.nanmin(avg_sig_pd, axis=2)

# Get NMSE at eta = 0.1
eta01_pow_nmse = avg_pow_nmse[-1, -3, :]
eta01_sig_nmse = avg_sig_nmse[-1, -3, :]

#####

fig, ax = plt.subplots()

linestyles = ['-', '--', '-.']
colors = ['tab:blue', 'tab:orange', 'tab:green']
markers = ['v', 'd']

a1, = ax.plot(1 - eta_range, worst_pow_pd[0, :], linewidth=0.0, linestyle=linestyles[0], color='black', marker=markers[0], label=r'power-based HRIS')
a2, = ax.plot(1 - eta_range, worst_sig_pd[0, :], linewidth=0.0, linestyle=linestyles[0], color='black', marker=markers[1], label=r'signal-based HRIS')

a3, = ax.plot(1 - eta_range, worst_pow_pd[0, :], linewidth=1.5, linestyle=linestyles[0], color=colors[0], label=r'$P_{\rm FA}=10^{-1}$')
a4, = ax.plot(1 - eta_range, worst_pow_pd[1, :], linewidth=1.5, linestyle=linestyles[1], color=colors[1], label=r'$P_{\rm FA}=10^{-2}$')
a5, = ax.plot(1 - eta_range, worst_pow_pd[2, :], linewidth=1.5, linestyle=linestyles[2], color=colors[2], label=r'$P_{\rm FA}=10^{-3}$')

# Go through all possible values for probability of false alarm
for pp, proba_false_alarm in enumerate(proba_false_alarm_range):

        ax.plot(1 - eta_range, worst_pow_pd[pp, :], linewidth=1.5, linestyle=linestyles[pp], marker=markers[0], color=colors[pp])
        ax.plot(1 - eta_range, worst_sig_pd[pp, :], linewidth=1.5, linestyle=linestyles[pp], marker=markers[1], color=colors[pp])

ax.set_xlabel(r'fraction of power absorbed, $1-\eta$')
ax.set_ylabel(r'probability of detection, $P_{\rm D}$')

ax.set_xscale('log')

ax.set_ylim([0.55, 1.05])

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.legend(framealpha=0.5)

plt.tight_layout()

plt.show()