import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

from scipy.interpolate import CubicSpline

import tikzplotlib

rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Load data
data = np.load('data/mmimo_K4.npz')

gen_data = np.load('data/figure7_gen-ris_K4_N32.npz')
pow_data = np.load('data/figure7_pow-ris_K4_N32.npz')
sig_data = np.load('data/figure7_sig-ris_K4_N32.npz')

# Extract relevant information
relative_probe_time = pow_data['relative_probe_time']

se = data['avg_se']

gen_se = gen_data['gen_avg_se']
pow_se = pow_data['pow_avg_se']
sig_se = sig_data['sig_avg_se']

# Get statistics
avg_se = se.mean(axis=(1, 2, 3))

avg_gen_se = gen_se.mean(axis=(2, 3, 4))
avg_pow_se = pow_se.mean(axis=(2, 3, 4))
avg_sig_se = sig_se.mean(axis=(2, 3, 4))

uq_se = np.percentile(se, q=75, axis=(1, 2, 3))

uq_gen_se = np.percentile(gen_se, q=75, axis=(2, 3, 4))
uq_pow_se = np.percentile(pow_se, q=75, axis=(2, 3, 4))
uq_sig_se = np.percentile(sig_se, q=75, axis=(2, 3, 4))

lq_se = np.percentile(se, q=25, axis=(1, 2, 3))

lq_gen_se = np.percentile(gen_se, q=25, axis=(2, 3, 4))
lq_pow_se = np.percentile(pow_se, q=25, axis=(2, 3, 4))
lq_sig_se = np.percentile(sig_se, q=25, axis=(2, 3, 4))

# Smoothness
st_relative_probe_time = np.linspace(relative_probe_time.min(), relative_probe_time.max(), 1001)

cs = CubicSpline(relative_probe_time[:-1], avg_gen_se[1, :-1])
st_gen_avg_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, avg_pow_se[1])
st_pow_avg_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, avg_sig_se[1])
st_sig_avg_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time[:-1], uq_gen_se[1, :-1])
st_gen_uq_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, uq_pow_se[1])
st_pow_uq_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, uq_sig_se[1])
st_sig_uq_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time[:-1], lq_gen_se[1, :-1])
st_gen_lq_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, lq_pow_se[1])
st_pow_lq_se = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, lq_sig_se[1])
st_sig_lq_se = cs(st_relative_probe_time)

linestyles = ['-', '--', '-.']
colors = ['tab:blue', 'tab:orange', 'tab:green']

# Plot figure for ZF
fig, ax = plt.subplots()

ax.plot(relative_probe_time, avg_se[1] * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='mMIMO baseline')
ax.fill_between(relative_probe_time, lq_se[1] * np.ones_like(relative_probe_time), uq_se[1] * np.ones_like(relative_probe_time), linewidth=0.0, color='black', alpha=0.1)

ax.plot(st_relative_probe_time, st_gen_avg_se, linewidth=1.5, linestyle='--', color=colors[0], label='genie baseline')
ax.fill_between(st_relative_probe_time, st_gen_lq_se, st_gen_uq_se, linewidth=0.0, color=colors[0], alpha=0.1)

ax.plot(st_relative_probe_time, st_pow_avg_se, linewidth=1.5, linestyle='-.', color=colors[1], label='power-based HRIS')
ax.fill_between(st_relative_probe_time, st_pow_lq_se, st_pow_uq_se, linewidth=0.0, color=colors[1], alpha=0.1)

ax.plot(st_relative_probe_time, st_sig_avg_se, linewidth=1.5, linestyle=':', color=colors[2], label='signal-based HRIS')
ax.fill_between(st_relative_probe_time, st_sig_lq_se, st_sig_uq_se, linewidth=0.0, color=colors[2], alpha=0.1)

ax.set_xlabel(r'relative probe duration, $\frac{C}{L}$')
ax.set_ylabel(r'$\underline{\mathrm{SE}}^{\rm UL}$ [bits/s/Hz]')

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

plt.tight_layout()

tikzplotlib.save("tikz/figure7c.tex")

plt.show()