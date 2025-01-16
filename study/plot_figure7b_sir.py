import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

from scipy.interpolate import CubicSpline

#import tikzplotlib

rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Load data
data = np.load('data/figure7_mmimo_K4.npz')

gen_data = np.load('data/figure7_gen-ris_K4_N32.npz')
pow_data = np.load('data/figure7_pow-ris_K4_N32.npz')
sig_data = np.load('data/figure7_sig-ris_K4_N32.npz')
breakpoint()
# Extract relevant infomation
relative_probe_time = pow_data['relative_probe_time']

avg_sir = data['avg_sir']

gen_avg_sir = gen_data['gen_avg_sir']
pow_avg_sir = pow_data['pow_avg_sir']
sig_avg_sir = sig_data['sig_avg_sir']

# Process
avg_sir = 10 * np.log10(avg_sir.mean(axis=(1, 2, 3)))

gen_avg_sir = 10 * np.log10(gen_avg_sir.mean(axis=(2, 3, 4)))
pow_avg_sir = 10 * np.log10(pow_avg_sir.mean(axis=(2, 3, 4)))
sig_avg_sir = 10 * np.log10(sig_avg_sir.mean(axis=(2, 3, 4)))

# Smoothness
st_relative_probe_time = np.linspace(relative_probe_time.min(), relative_probe_time.max(), 1001)

cs = CubicSpline(relative_probe_time[:-1], gen_avg_sir[1, :-1])
st_gen_avg_sir = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, pow_avg_sir[1])
st_pow_avg_sir = cs(st_relative_probe_time)

cs = CubicSpline(relative_probe_time, sig_avg_sir[1])
st_sig_avg_sir = cs(st_relative_probe_time)

fig, ax = plt.subplots()

ax.plot(relative_probe_time, avg_sir[1] * np.ones_like(relative_probe_time), linewidth=1.5, linestyle='-', color='black', label='MMIMO')

ax.plot(st_relative_probe_time, st_gen_avg_sir, linewidth=1.5, linestyle='--', label='Genie baseline')
ax.plot(st_relative_probe_time, st_pow_avg_sir, linewidth=1.5, linestyle='-.', label='Power-based HRIS')
ax.plot(st_relative_probe_time, st_sig_avg_sir, linewidth=1.5, linestyle=':', label='Signal-based HRIS')

ax.set_xlabel(r'relative probe duration, $\frac{C}{L}$')
ax.set_ylabel(r'$\underline{\mathrm{SIR}}_k$ [dB]')

ax.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

#tikzplotlib.save("tikz/figure7b.tex")

plt.tight_layout()

plt.show()