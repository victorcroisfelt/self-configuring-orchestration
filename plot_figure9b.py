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

ctl_data = np.load('data/figure7_control-ris_K4_N32.npz')

# Number of pilot subbblocks for probe
n_pilot_subblocks = int(64 // 4)
n_pilot_subblocks_probe_range = np.arange(1, n_pilot_subblocks + 1)
X_ = n_pilot_subblocks_probe_range / n_pilot_subblocks 

# Extract relevant information
se = data['avg_se_pos']

gen_se = gen_data['gen_avg_se_pos']#[:, :8]
pow_se = pow_data['pow_avg_se_pos']#[:, :8]
sig_se = sig_data['sig_avg_se_pos']#[:, :8]

ctl_se = ctl_data['avg_se']

#breakpoint()
# Get statistics
avg_se = se.mean(axis=(1, 2, 3, 4))

avg_gen_se = np.nanmean(gen_se, axis=(2, 3, 4, 5))
avg_pow_se = np.nanmean(pow_se, axis=(2, 3, 4, 5))
avg_sig_se = np.nanmean(sig_se, axis=(2, 3, 4, 5))

avg_ctl_se = np.nanmean(ctl_se, axis=(1, 2, 3, 4))

# uq_se = np.nanpercentile(se, q=75, axis=(1, 2, 3, 4))
#
# uq_gen_se = np.nanpercentile(gen_se, q=75, axis=(2, 3, 4, 5))
# uq_pow_se = np.nanpercentile(pow_se, q=75, axis=(2, 3, 4, 5))
# uq_sig_se = np.nanpercentile(sig_se, q=75, axis=(2, 3, 4, 5))
#
# lq_se = np.nanpercentile(se, q=25, axis=(1, 2, 3, 4))
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
fig, ax = plt.subplots(figsize=(2*1.5, 2.0))

ax.plot(X_, avg_se[0] * np.ones_like(X_), linewidth=1.5, linestyle='-', color='black', label='mMIMO baseline')
#ax.fill_between(X_, lq_se[0] * np.ones_like(X_), uq_se[0] * np.ones_like(X_), linewidth=0.0, color='black', alpha=0.1)

ax.plot(X_, avg_gen_se[0], linewidth=1.5, linestyle='--', markersize=4, marker=markers[0], color=colors[0], label='Genie baseline')
#ax.fill_between(X_, lq_gen_se[0], uq_gen_se[0], linewidth=0.0, color=colors[0], alpha=0.1)

ax.plot(X_, avg_pow_se[0], linewidth=1.5, linestyle='-.', markersize=4, marker=markers[1], color=colors[1], label='PD-enabled HRIS')
#ax.fill_between(X_, lq_pow_se[0], uq_pow_se[0], linewidth=0.0, color=colors[1], alpha=0.1)

ax.plot(X_, avg_sig_se[0], linewidth=1.5, linestyle=':', markersize=4, marker=markers[2], color=colors[2], label='DSP-enabled HRIS')
#ax.fill_between(X_, lq_se[0], uq_sig_se[0], linewidth=0.0, color=colors[2], alpha=0.1)

ax.plot(X_, (128 - 72)/128 * avg_ctl_se[0] * np.ones_like(X_), linewidth=1.5, label='Controled RIS')

controled_ris_best = (128 - 40)/128 * avg_ctl_se[0]
controled_ris_worse = (128 - 72)/128 * avg_ctl_se[0]

ax.set_xlabel(r'Relative probe duration within CHEST, $\varpi$', fontsize=10)
ax.set_ylabel(r'$\underline{\mathrm{SE}}^{\rm UL}$ [bits/s/Hz]', fontsize=10)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

ax.grid(True, linestyle='--', linewidth=0.25, color='gray', alpha=0.5)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

plt.tight_layout()

ax.legend(fontsize=8)

#--------------------
# Save data for TiKz plotting
#--------------------

# Informed BS
Y_ = avg_gen_se[0]
informed = np.vstack((X_, Y_)).T
np.savetxt('txts/fig9_b_informed.txt', informed, fmt='%.4f', comments='')

# PD-enabled HRIS 
Y_ = avg_pow_se[0]
pd = np.vstack((X_, Y_)).T
np.savetxt('txts/fig9_b_pd.txt', pd, fmt='%.4f', comments='')

# DSP-enabled HRIS
Y_ = avg_sig_se[0]
dsp = np.vstack((X_, Y_)).T
np.savetxt('txts/fig9_b_dsp.txt', dsp, fmt='%.4f', comments='')