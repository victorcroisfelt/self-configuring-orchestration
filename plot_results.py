import scipy.io

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc

# LaTeX type definitions
rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

########################################
# Load data
########################################
mat = scipy.io.loadmat('data/RESULTS_M64_N32_K16_L64_250x250.mat')

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

# Flatten results in terms of setup and channel realizations
th_detected_ue_pow = th_detected_ue_pow.reshape(*th_detected_ue_pow.shape[:2], -1)
th_detected_ue_sig = th_detected_ue_sig.reshape(*th_detected_ue_sig.shape[:2], -1)

hat_detected_ue_pow = hat_detected_ue_pow.reshape(*hat_detected_ue_pow.shape[:2], -1)
hat_detected_ue_sig = hat_detected_ue_sig.reshape(*hat_detected_ue_sig.shape[:2], -1)

mse_cha_est_pow = mse_cha_est_pow.reshape(*mse_cha_est_pow.shape[:2], -1)
mse_cha_est_sig = mse_cha_est_sig.reshape(*mse_cha_est_sig.shape[:2], -1)

########################################
# Plot
########################################
fig, ax = plt.subplots()

linestyles = ['-', '--', '-.']

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(false_alarm_prob_vec):

 
	dumb = ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

	ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color='tab:blue')
	ax.plot(C_vec/L, hat_detected_ue_sig[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color='tab:orange')

	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

ax.set_xlabel(r'Relative duration $(C/L)$')
ax.set_ylabel(r'Probability of detection $P_{\rm D}$')

ax.legend()

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

dumb.pop().remove()

plt.show()


fig, ax = plt.subplots()

linestyles = ['-', '--', '-.']

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(false_alarm_prob_vec):

	dumb = ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

	ax.plot(C_vec/L, mse_cha_est_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color='tab:blue')
	ax.plot(C_vec/L, mse_cha_est_sig[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color='tab:orange')

	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

ax.set_xlabel(r'Relative duration $(C/L)$')
ax.set_ylabel('MSE')

ax.set_yscale('log')

ax.legend()

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

dumb.pop().remove()

plt.show()


# figure
# hold on
# p11 = plot(C_vec./L,MSE_cha_est_pow(:,1));
# p12 = plot(C_vec./L,MSE_cha_est_sig(:,1));
# % p21 = plot(C_vec./L,MSE_cha_est_pow(:,2));
# % p22 = plot(C_vec./L,MSE_cha_est_sig(:,2));
# % p31 = plot(C_vec./L,MSE_cha_est_pow(:,3));
# % p32 = plot(C_vec./L,MSE_cha_est_sig(:,3));
# set(gca, 'YScale', 'log')


# ylabel('Channel estiation MSE', 'Interpreter', 'latex')
# xlabel('$C/L$', 'Interpreter', 'latex')

# legend('power-based', 'signal-based', 'Interpreter', 'latex')


# figure
# hold on
# p11 = plot(C_vec./L,avg_hat_detected_ue_pow(:,1));
# p12 = plot(C_vec./L,avg_hat_detected_ue_sig(:,1));
# p21 = plot(C_vec./L,avg_hat_detected_ue_pow(:,2));
# p22 = plot(C_vec./L,avg_hat_detected_ue_sig(:,2));
# p31 = plot(C_vec./L,avg_hat_detected_ue_pow(:,3));
# p32 = plot(C_vec./L,avg_hat_detected_ue_sig(:,3));


# ylabel('Detecion probability', 'Interpreter', 'latex')
# xlabel('$C/L$', 'Interpreter', 'latex')

# legend(['power-based ', num2str(false_alarm_prob_vec(1))], ['signal-based'], 'Interpreter', 'latex')
