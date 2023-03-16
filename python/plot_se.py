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

avg_se = data['avg_se']
pow_avg_se = pow_data['pow_avg_se']
sig_avg_se = sig_data['sig_avg_se']

fig, ax = plt.subplots()

ax.plot(relative_probe_time, avg_se.mean() * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='MMIMO')
ax.plot(relative_probe_time, pow_avg_se.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle='--', label='Power-based HRIS')
ax.plot(relative_probe_time, sig_avg_se.mean(axis=(1, 2, 3)), linewidth=1.5, linestyle=':', label='Signal-based HRIS')

ax.set_xlabel(r'relative probe duration $\frac{C}{L}$')
ax.set_ylabel(r'avg. SE per UE [bits/s/Hz]')

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

ax.legend()

plt.tight_layout()

plt.show()