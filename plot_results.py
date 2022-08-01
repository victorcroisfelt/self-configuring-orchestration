import scipy.io

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc

# LaTeX type definitions
rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern'], 'size': 14})
#rc('text', usetex=True)

########################################
# Load data
########################################
mat = scipy.io.loadmat('data/RESULTS_M64_N16_K16_L128_250x250_mod.mat')

# Loading parameters 
L = np.squeeze(mat['L'])
C_vec = np.squeeze(mat['C_vec'])
false_alarm_prob_vec = np.squeeze(mat['false_alarm_prob_vec'])

# Loading results 
th_detected_ue_pow = mat['th_detected_ue_pow']
th_detected_ue_sig = mat['th_detected_ue_sig']

hat_detected_ue_pow = mat['hat_detected_ue_pow']
hat_detected_ue_sig = mat['hat_detected_ue_sig']

mse_cha_est_pow = mat['MSE_cha_est_pow']
mse_cha_est_sig = mat['MSE_cha_est_sig']

sinr_est_pow = mat['SINR_est_sig']
sinr_est_sig = mat['SINR_est_pow']

se_est_pow = mat['SE_est_sig']
se_est_sig = mat['SE_est_pow']

# Flatten results in terms of setup and channel realizations
th_detected_ue_pow = th_detected_ue_pow.reshape(*th_detected_ue_pow.shape[:2], -1)
th_detected_ue_sig = th_detected_ue_sig.reshape(*th_detected_ue_sig.shape[:2], -1)

hat_detected_ue_pow = hat_detected_ue_pow.reshape(*hat_detected_ue_pow.shape[:2], -1)
hat_detected_ue_sig = hat_detected_ue_sig.reshape(*hat_detected_ue_sig.shape[:2], -1)

mse_cha_est_pow = mse_cha_est_pow.reshape(*mse_cha_est_pow.shape[:2], -1)
mse_cha_est_sig = mse_cha_est_sig.reshape(*mse_cha_est_sig.shape[:2], -1)

sinr_est_pow = sinr_est_pow.reshape(*sinr_est_pow.shape[:3], -1).mean(axis=0)
sinr_est_sig = sinr_est_sig.reshape(*sinr_est_sig.shape[:3], -1).mean(axis=0)

se_est_pow = se_est_pow.reshape(*se_est_pow.shape[:3], -1).mean(axis=0)
se_est_sig = se_est_sig.reshape(*se_est_sig.shape[:3], -1).mean(axis=0)

########################################
# Plot
########################################
fig_size = (6.5, 3.5)

fig, ax = plt.subplots(figsize=fig_size)

#linestyles = ['-', '--', '-.']
markers = ['.', '+', 'x']
dumb_list = []

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(false_alarm_prob_vec):

 
	dumb = ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=0.0, marker=markers[fa], color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

	ax.plot(C_vec/L, hat_detected_ue_sig[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle='-', marker=markers[fa], color='tab:blue')
	ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle='--', marker=markers[fa], color='tab:orange')
	

	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	dumb_list.append(dumb)

dumb = ax.plot(C_vec/L, hat_detected_ue_sig[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle='-', color='tab:blue', label='Signal-based')
dumb_list.append(dumb)
dumb = ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle='--', color='tab:orange', label='Power-based')
dumb_list.append(dumb)

ax.set_xlabel(r'Relative duration $(C/L)$')
ax.set_ylabel(r'Probability of detection $P_{\rm D}$')

ax.legend(fontsize='x-small', framealpha=0.5)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

[(dumb.pop(0)).remove() for dumb in dumb_list]

plt.tight_layout()

plt.savefig('figures/detection.pdf', dpi='figure', format='pdf', transparent='True')

#####


fig, ax = plt.subplots(figsize=fig_size)

dumb_list = []

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(false_alarm_prob_vec):

	dumb = ax.plot(C_vec/L, mse_cha_est_pow[:, fa, :].mean(axis=-1), linewidth=0.0, marker=markers[fa], color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

	ax.plot(C_vec/L, mse_cha_est_sig[:, fa, :].mean(axis=-1), linewidth=1.5, marker=markers[fa], linestyle='-', color='tab:blue')
	ax.plot(C_vec/L, mse_cha_est_pow[:, fa, :].mean(axis=-1), linewidth=1.5, marker=markers[fa], linestyle='--', color='tab:orange')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	dumb_list.append(dumb)

ax.set_xlabel(r'Relative duration $(C/L)$')
ax.set_ylabel('MSE')

ax.set_yscale('log')

#ax.legend(fontsize='x-small', framealpha=0.5)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

[(dumb.pop(0)).remove() for dumb in dumb_list]

plt.tight_layout()

plt.savefig('figures/mse.pdf', dpi='figure', format='pdf', transparent='True')


#####


fig, ax = plt.subplots(figsize=fig_size)

dumb_list = []

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(false_alarm_prob_vec):

	dumb = ax.plot(C_vec/L, 10 * np.log10(sinr_est_pow[:, fa, :].mean(axis=-1)), linewidth=0.0, marker=markers[fa], color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

	ax.plot(C_vec/L, 10 * np.log10(sinr_est_sig[:, fa, :].mean(axis=-1)), linewidth=1.5, marker=markers[fa], linestyle='-', color='tab:blue')
	ax.plot(C_vec/L, 10 * np.log10(sinr_est_pow[:, fa, :].mean(axis=-1)), linewidth=1.5, marker=markers[fa], linestyle='--', color='tab:orange')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	dumb_list.append(dumb)

ax.set_xlabel(r'Relative duration $(C/L)$')
ax.set_ylabel('SINR [dB]')

#ax.legend(fontsize='x-small', framealpha=0.5)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

[(dumb.pop(0)).remove() for dumb in dumb_list]

plt.tight_layout()

plt.savefig('figures/sinr.pdf', dpi='figure', format='pdf', transparent='True')

#####


fig, ax = plt.subplots(figsize=fig_size)

dumb_list = []

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(false_alarm_prob_vec):

	dumb = ax.plot(C_vec/L, mse_cha_est_pow[:, fa, :].mean(axis=-1), linewidth=0.0, marker=markers[fa], color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

	ax.plot(C_vec/L, se_est_sig[:, fa, :].mean(axis=-1), linewidth=1.5, marker=markers[fa], linestyle='-', color='tab:blue')
	ax.plot(C_vec/L, se_est_pow[:, fa, :].mean(axis=-1), linewidth=1.5, marker=markers[fa], linestyle='--', color='tab:orange')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	dumb_list.append(dumb)

ax.set_xlabel(r'Relative duration $(C/L)$')
ax.set_ylabel('SE [bit/s/Hz]')

ax.legend(fontsize='x-small', framealpha=0.5)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

[(dumb.pop(0)).remove() for dumb in dumb_list]

plt.tight_layout()

#plt.show()


