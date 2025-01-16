import matplotlib.pyplot as plt
import numpy as np

# Get data
data = np.load('data/figure6a.npz')

pow_test = data['pow_test'].reshape((4, 16, -1))
pow_threshold = data['pow_threshold']

sig_test = data['sig_test'].reshape((4, 16, -1))
sig_threshold = data['sig_threshold']

# Evaluate how
avg_pow_test = pow_test.mean(axis=-1)
avg_sig_test = sig_test.mean(axis=-1)

#####
# Effect of eta
#####
sort_eta_pow = np.sort(pow_test, axis=-1)
sort_eta_sig = np.sort(sig_test, axis=-1)

y_axis = np.linspace(0, 1,40000)

fig, axes = plt.subplots(ncols=2)

for rr in range(4):
    axes[0].plot(sort_eta_pow[rr, 0], y_axis)
    axes[1].plot(sort_eta_sig[rr, 0], y_axis)

axes[0].plot([pow_threshold[0], pow_threshold[0]], [0, 1], color='black')
axes[1].plot([sig_threshold[0], sig_threshold[0]], [0, 1], color='black')

axes[0].set_xscale('log')
axes[1].set_xscale('log')

breakpoint()