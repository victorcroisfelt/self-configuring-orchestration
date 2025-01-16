import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

from scipy.interpolate import CubicSpline

#import tikzplotlib

# Configure matplotlib to use the Helvetica Neue font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 10

# Load data
data = np.load('data/figure7_mmimo_K4.npz')

gen_data = np.load('data/figure7_gen-ris_K4_N32.npz')
pow_data = np.load('data/figure7_pow-ris_K4_N32.npz')
sig_data = np.load('data/figure7_sig-ris_K4_N32.npz')

# Number of pilot subbblocks for probe
X_ = (16 - np.arange(1,  16 + 1)) / 16

gen_se = gen_data['gen_avg_se'][:, -2]
pow_se = pow_data['pow_avg_se_pos'][:, 8]
sig_se = sig_data['sig_avg_se_pos'][:, 8]


avg_gen_se = np.nanmean(gen_se, axis=(1, 2, 3, 4))
avg_pow_se = np.nanmean(pow_se, axis=(1, 2, 3, 4))
avg_sig_se = np.nanmean(sig_se, axis=(1, 2, 3, 4))

# Number of pilot subbblocks for probe
relative_control_duration = np.linspace(0, 0.1)

# Coherence interval length
tau_c = 128

# Calculate pre-log term
pre_log_term = (tau_c - 16 * 4 - relative_control_duration * tau_c) / tau_c
avg_gen_se = avg_gen_se[:, None] * pre_log_term[None, :]

breakpoint()

# uq_gen_se = np.nanpercentile(gen_se, q=75, axis=(2, 3, 4, 5))
# uq_pow_se = np.nanpercentile(pow_se, q=75, axis=(2, 3, 4, 5))
# uq_sig_se = np.nanpercentile(sig_se, q=75, axis=(2, 3, 4, 5))
#
# lq_gen_se = np.nanpercentile(gen_se, q=25, axis=(2, 3, 4, 5))
# lq_pow_se = np.nanpercentile(pow_se, q=25, axis=(2, 3, 4, 5))
# lq_sig_se = np.nanpercentile(sig_se, q=25, axis=(2, 3, 4, 5))

# Smoothness
# st_relative_probe_time = np.linspace(relative_probe_time.min(), relative_probe_time.max(), 1001)
#
# cs = CubicSpline(relative_probe_time[:-1], avg_gen_se[1, :-1])
# st_gen_avg_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time, avg_pow_se[1])
# st_pow_avg_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time, avg_sig_se[1])
# st_sig_avg_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time[:-1], uq_gen_se[1, :-1])
# st_gen_uq_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time, uq_pow_se[1])
# st_pow_uq_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time, uq_sig_se[1])
# st_sig_uq_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time[:-1], lq_gen_se[1, :-1])
# st_gen_lq_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time, lq_pow_se[1])
# st_pow_lq_se = cs(st_relative_probe_time)
#
# cs = CubicSpline(relative_probe_time, lq_sig_se[1])
# st_sig_lq_se = cs(st_relative_probe_time)

linestyles = ['-', '--', '-.']
colors = ['tab:blue', 'tab:orange', 'tab:green']
markers = ['*', 'p', 'd']

# Plot figure for ZF
fig, ax = plt.subplots(figsize=(2*1.5, 4))
ax.plot(relative_control_duration, avg_gen_se[1], linewidth=1.5, linestyle='--', markersize=4, marker=markers[0], color=colors[0], label='Controlled RIS')
ax.plot(relative_control_duration, avg_pow_se[1] * np.ones_like(relative_control_duration), linewidth=1.5, linestyle='-.', markersize=4, marker=markers[1], color=colors[1], label='Power-based HRIS')
ax.plot(relative_control_duration, avg_sig_se[1] * np.ones_like(relative_control_duration), linewidth=1.5, linestyle=':', markersize=4, marker=markers[2], color=colors[2], label='Signal-based HRIS')
#ax.fill_between(X_, lq_se[0], uq_sig_se[0], linewidth=0.0, color=colors[2], alpha=0.1)

ax.set_xlabel(r'Relative control duration', fontsize=10)
ax.set_ylabel(r'$\underline{\mathrm{SE}}^{\rm UL}$ [bits/s/Hz]', fontsize=10)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

ax.grid(True, linestyle='--', linewidth=0.25, color='gray', alpha=0.5)

ax.set_xticks([0, 0.1])

plt.tight_layout()

#tikzplotlib.save("tikz/figure7c.tex")
plt.savefig('figs/figure10.pdf', bbox_inches='tight')

plt.show()