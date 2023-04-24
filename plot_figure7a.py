import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline

from matplotlib import rc

import tikzplotlib

rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Load data
data = np.load('data/figure7_mmimo_K4.npz')

gen_data = np.load('data/figure7_gen-ris_K4_N32.npz')
pow_data = np.load('data/figure7_pow-ris_K4_N32.npz')
sig_data = np.load('data/figure7_sig-ris_K4_N32.npz')

# Extract relevant information
relative_probe_time = pow_data['relative_probe_time']

nmse = data['avg_nmse']

gen_nmse = gen_data['gen_avg_nmse']
pow_nmse = pow_data['pow_avg_nmse']
sig_nmse = sig_data['sig_avg_nmse']

# Get statistics
avg_nmse = nmse.mean()

avg_gen_nmse = gen_nmse.mean(axis=(1, 2, 3))
avg_pow_nmse = pow_nmse.mean(axis=(1, 2, 3))
avg_sig_nmse = sig_nmse.mean(axis=(1, 2, 3))

# Smoothness
st_relative_probe_time = np.linspace(relative_probe_time.min(), relative_probe_time.max(), 1001)

cs = CubicSpline(relative_probe_time[:-1], avg_gen_nmse[:-1])
st_avg_gen_nmse = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, avg_pow_nmse)
st_avg_pow_nmse = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, avg_sig_nmse)
st_avg_sig_nmse = cs(st_relative_probe_time)

# Plot figure for MR
fig, ax = plt.subplots()

#linestyles = ['-', '--', '-.']
colors = ['tab:blue', 'tab:orange', 'tab:green']
#markers = ['*', 'v', 'd']

ax.plot(relative_probe_time, avg_nmse * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='mMIMO baseline')
ax.plot(st_relative_probe_time, st_avg_gen_nmse, linewidth=1.5, linestyle='--', color=colors[0], label='genie baseline')
ax.plot(st_relative_probe_time, st_avg_pow_nmse, linewidth=1.5, linestyle='-.', color=colors[1], label='power-based HRIS')
ax.plot(st_relative_probe_time, st_avg_sig_nmse, linewidth=1.5, linestyle=':', color=colors[2], label='signal-based HRIS')

ax.set_xlabel(r'relative probe duration, $\frac{C}{L}$')
ax.set_ylabel(r'avg. $\mathrm{NMSE}_{\mathrm{CHEST}, k}$ in (45)')

ax.set_yscale('log')

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

ax.legend(framealpha=0.5)

plt.tight_layout()

tikzplotlib.save("tikz/figure7a.tex")

plt.show()