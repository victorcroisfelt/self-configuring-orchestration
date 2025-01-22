import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline

from matplotlib import rc

#import tikzplotlib

# Configure matplotlib to use the Helvetica Neue font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 10

# Load data
data = np.load('data/figure9_mmimo_K4.npz')

gen_data = np.load('data/figure9_gen-ris_K4_N32.npz')
pow_data = np.load('data/figure9_pow-ris_K4_N32.npz')
sig_data = np.load('data/figure9_sig-ris_K4_N32.npz')

# Extract relevant information
nmse = data['avg_nmse']

gen_nmse = gen_data['gen_avg_nmse']
pow_nmse = pow_data['pow_avg_nmse']
sig_nmse = sig_data['sig_avg_nmse']

gen_nmse = gen_nmse[:16]
pow_nmse = pow_nmse[:16]
sig_nmse = sig_nmse[:16]

# Get statistics
avg_nmse = nmse.mean()

avg_gen_nmse = np.nanmean(gen_nmse, axis=(1, 2, 3, 4))
avg_pow_nmse = np.nanmean(pow_nmse, axis=(1, 2, 3, 4))
avg_sig_nmse = np.nanmean(sig_nmse, axis=(1, 2, 3, 4))

# Number of pilot subbblocks for probe
X_ = np.arange(1,  16 + 1) / 16

#linestyles = ['-', '--', '-.']
colors = ['tab:blue', 'tab:orange', 'tab:green']
markers = ['*', 'p', 'd']

# Plot figure for MR
fig, ax = plt.subplots(figsize=(2*1.5, 2.0))

ax.plot(X_, avg_nmse * np.ones_like(X_), linewidth=1.0, linestyle='-', color='black', label='mMIMO baseline')

ax.plot(X_, avg_gen_nmse, linewidth=1.0, linestyle='--', color=colors[0], markersize=4, marker=markers[0], label='Genie BS baseline')
ax.plot(X_, avg_pow_nmse, linewidth=1.0, linestyle='-.', color=colors[1], markersize=4, marker=markers[1], label='Power-based HRIS')
ax.plot(X_, avg_sig_nmse, linewidth=1.0, linestyle=':', color=colors[2], markersize=4, marker=markers[2], label='Signal-based HRIS')

ax.set_xlabel(r'Rel. reflec. duration during CHEST, $\varpi$', fontsize=10)
ax.set_ylabel(r'avg. $\mathrm{NMSE}_{\mathrm{CHEST}, k}$', fontsize=10)

#ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
ax.set_yticks([10e-3, 10e-4, 10e-5])

ax.set_yscale('log')

# ax.set_xscale('log')
#
#ax.set_xticklabels(ax.get_xticks(), fontsize=8)
#ax.set_yticklabels(ax.get_yticks(), fontsize=8)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

ax.grid(True, linestyle='--', linewidth=0.25, color='gray', alpha=0.5)

legend = ax.legend(fontsize=8)
#legend.get_frame().set_facecolor('lightgray')
legend.get_frame().set_linewidth(0.5)
plt.tight_layout()

#tikzplotlib.save("tikz/figure7a.tex")
#plt.savefig('figs/figure9a.pdf', bbox_inches='tight')

plt.show()