import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Load data
data = np.load('data/mmimo_K4.npz')
pow_data = np.load('data/pow-ris_K4.npz')
sig_data = np.load('data/sig-ris_K4.npz')

# Extract relevant infomation
relative_probe_time = pow_data['relative_probe_time']

avg_num = 10 * np.log10(data['avg_num'])
pow_avg_num = 10 * np.log10(pow_data['pow_avg_num'])
sig_avg_num = 10 * np.log10(sig_data['sig_avg_num'])

avg_den1 = 10 * np.log10(data['avg_den1'])
pow_avg_den1 = 10 * np.log10(pow_data['pow_avg_den1'])
sig_avg_den1 = 10 * np.log10(sig_data['sig_avg_den1'])

avg_den2 = 10 * np.log10(data['avg_den2'])
pow_avg_den2 = 10 * np.log10(pow_data['pow_avg_den2'])
sig_avg_den2 = 10 * np.log10(sig_data['sig_avg_den2'])

fig, axes = plt.subplots(ncols=3)

axes[0].plot(relative_probe_time, avg_num.mean() * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='MMIMO')
axes[0].plot(relative_probe_time, pow_avg_num.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle='--', label='Power-based HRIS')
axes[0].plot(relative_probe_time, sig_avg_num.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle=':', label='Signal-based HRIS')

axes[0].set_xlabel(r'relative probe duration $\frac{C}{L}$')
axes[0].set_ylabel(r'num. [dB]')

axes[0].grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

axes[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

axes[0].legend(fontsize='x-small')


axes[1].plot(relative_probe_time, avg_den1.mean() * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='MMIMO')
axes[1].plot(relative_probe_time, pow_avg_den1.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle='--', label='Power-based HRIS')
axes[1].plot(relative_probe_time, sig_avg_den1.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle=':', label='Signal-based HRIS')

axes[1].set_xlabel(r'relative probe duration $\frac{C}{L}$')
axes[1].set_ylabel(r'denom. 1 (interf.) [dB]')

axes[1].grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

axes[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])


axes[2].plot(relative_probe_time, avg_den2.mean() * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='MMIMO')
axes[2].plot(relative_probe_time, pow_avg_den2.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle='--', label='Power-based HRIS')
axes[2].plot(relative_probe_time, sig_avg_den2.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle=':', label='Signal-based HRIS')

axes[2].set_xlabel(r'relative probe duration $\frac{C}{L}$')
axes[2].set_ylabel(r'denom. 2 (noise) [dB]')

axes[2].grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

axes[2].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

plt.tight_layout()

plt.show()