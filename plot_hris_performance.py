import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

#import tikzplotlib

# # Find the path to the Helvetica Neue font file
# font_path = fm.findfont(fm.FontProperties(family='Helvetica Neue'))

# Configure matplotlib to use the Helvetica Neue font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 10

# Get data
data = np.load('data/figure6a.npz')

pow_detection = data['pow_detection'].reshape((2, 4, 16, -1))
sig_detection = data['sig_detection'].reshape((2, 4, 16, -1))

avg_pow_detection = pow_detection.mean(axis=-1)
avg_sig_detection = sig_detection.mean(axis=-1)

# Number of pilot subbblocks for probe
X_ = (16 - np.arange(1,  16 + 1)) / 16

# Range of HRIS reflection parameter
Y_ = 1 - np.array([0.9, 0.99, 0.999, 0.9999])

# Customization
linestyles = ['-', '--', '-.']
colors = ['tab:orange', 'tab:green']
markers = ['p', 'd']

fig, ax = plt.subplots(figsize=(2*1.5, 2.0))

a1, = ax.plot(Y_, avg_pow_detection[0, :, 8], linewidth=0.0, linestyle=linestyles[0], markersize=4, color='black', marker=markers[0], label=r'Power-based HRIS')
a2, = ax.plot(Y_, avg_sig_detection[0, :, 8], linewidth=0.0, linestyle=linestyles[0], markersize=4, color='black', marker=markers[1], label=r'Signal-based HRIS')

a3, = ax.plot(Y_, avg_pow_detection[0, :, 8], linewidth=1.0, linestyle=linestyles[0], markersize=4, color='black', label=r'$P_{\rm FA}=10^{-2}$')
a4, = ax.plot(Y_, avg_pow_detection[1, :, 8], linewidth=1.0, linestyle=linestyles[1], markersize=4, color='black', label=r'$P_{\rm FA}=10^{-3}$')
#a5, = ax.plot(Y_, avg_pow_detection[2, :, -1], linewidth=1.5, linestyle=linestyles[2], color=colors[2], label=r'$P_{\rm FA}=10^{-3}$')

ax.plot(Y_, avg_pow_detection[0, :, 8], linewidth=1.0, linestyle=linestyles[0], markersize=4, color=colors[0], marker=markers[0])
ax.plot(Y_, avg_pow_detection[1, :, 8], linewidth=1.0, linestyle=linestyles[1], markersize=4, color=colors[1], marker=markers[0])
#ax.plot(Y_, avg_pow_detection[2, :, -1], linewidth=1.5, linestyle=linestyles[2], color=colors[2], marker=markers[0])

ax.plot(Y_, avg_sig_detection[0, :, 8], linewidth=1.0, linestyle=linestyles[0], markersize=4, color=colors[0], marker=markers[1])
ax.plot(Y_, avg_sig_detection[1, :, 8], linewidth=1.0, linestyle=linestyles[1], markersize=4, color=colors[1], marker=markers[1])
#ax.plot(Y_, avg_sig_detection[2, :, -1], linewidth=1.5, linestyle=linestyles[2], color=colors[2], marker=markers[1])

ax.set_xlabel(r'Frac. of power abs. by the HRIS, $1-\eta$', fontsize=10)
ax.set_ylabel(r'Proba. of detection, $P_{\rm D}$', fontsize=10)

#ax.set_ylim([0.7, 1.005])

ax.set_xticks([10e-4, 10e-3, 10e-2, 10e-1])
ax.set_yticks([0.7, 0.8, 0.9, 1])

ax.set_xscale('log')

# ax.set_xticklabels(ax.get_xticks(), fontsize=8)
# ax.set_yticklabels(ax.get_yticks(), fontsize=8)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

ax.grid(True, linestyle='--', linewidth=0.25, color='gray', alpha=0.5)

legend = ax.legend(fontsize=8)
#legend.get_frame().set_facecolor('lightgray')
legend.get_frame().set_linewidth(0.5)

a1.remove()
a2.remove()
a3.remove()
a4.remove()
#a5.remove()

plt.tight_layout()

tikzplotlib.save("tikz/figure7a.tex")
plt.savefig('figs/figure7a.pdf', bbox_inches='tight')

#plt.show()


#####


fig, ax = plt.subplots(figsize=(2*1.5, 2.0))

ax.plot(X_, avg_pow_detection[0, -1, :], linewidth=1.0, markersize=4, linestyle=linestyles[0], color=colors[0], marker=markers[0])
ax.plot(X_, avg_pow_detection[1, -1, :], linewidth=1.0, markersize=4, linestyle=linestyles[1], color=colors[1], marker=markers[0])
#ax.plot(X_, avg_pow_detection[2, 1, :], linewidth=1.5, linestyle=linestyles[2], color=colors[2], marker=markers[0])

ax.plot(X_, avg_sig_detection[0, -1, :], linewidth=1.0, markersize=4, linestyle=linestyles[0], color=colors[0], marker=markers[1])
ax.plot(X_, avg_sig_detection[1, -1, :], linewidth=1.0, markersize=4, linestyle=linestyles[1], color=colors[1], marker=markers[1])
#ax.plot(X_, avg_sig_detection[2, 1, :], linewidth=1.5, linestyle=linestyles[2], color=colors[2], marker=markers[1])

ax.set_xlabel(r'Rel. reflec. duration during CHEST, $\varpi$', fontsize=10)
ax.set_ylabel(r'Proba. of detection, $P_{\rm D}$', fontsize=10)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1])

# ax.set_xticklabels(ax.get_xticks(), fontsize=8)
# ax.set_yticklabels(ax.get_yticks(), fontsize=8)

ax.grid(True, linestyle='--', linewidth=0.25, color='gray', alpha=0.5)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

#ax.legend(framealpha=0.5)
plt.tight_layout()

tikzplotlib.save("tikz/figure7b.tex")
plt.savefig('figs/figure7b.pdf', bbox_inches='tight')

plt.show()