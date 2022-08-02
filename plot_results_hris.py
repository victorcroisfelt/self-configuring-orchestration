import scipy.io

import numpy as np

import matplotlib.pyplot as plt

# LaTeX type definitions
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amssymb,amsmath,amsfonts,amsthm,mathtools,cuted,bbold} \usepackage[cmintegrals]{newtxmath}')

########################################
# Load data
########################################
filename = 'data/RESULTS_M64_N16_K16_L128_100x100.mat'

mat = scipy.io.loadmat(filename)

# Loading parameters 
L = np.squeeze(mat['L'])
C_vec = np.squeeze(mat['C_vec'])
false_alarm_prob_vec = np.squeeze(mat['false_alarm_prob_vec'])

# Loading results 
th_detected_ue_pow = mat['th_detected_ue_pow']
th_detected_ue_sig = mat['th_detected_ue_sig']

hat_detected_ue_pow = mat['hat_detected_ue_pow']
hat_detected_ue_sig = mat['hat_detected_ue_sig']

distance_pow = mat['distance_pow']
distance_sig = mat['distance_sig']

# Flatten results in terms of setup and channel realizations
th_detected_ue_pow = th_detected_ue_pow.reshape(*th_detected_ue_pow.shape[:2], -1)
th_detected_ue_sig = th_detected_ue_sig.reshape(*th_detected_ue_sig.shape[:2], -1)

hat_detected_ue_pow = hat_detected_ue_pow.reshape(*hat_detected_ue_pow.shape[:2], -1)
hat_detected_ue_sig = hat_detected_ue_sig.reshape(*hat_detected_ue_sig.shape[:2], -1)

distance_pow = distance_pow.reshape(*distance_pow.shape[:2], -1)
distance_sig = distance_sig.reshape(*distance_sig.shape[:2], -1)

########################################
# Plot
########################################
fig_size = (6.5/2, 3.5)
#xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
xticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

fig, ax = plt.subplots(figsize=fig_size)

linestyles = ['-', '--', '-.']
labels = [r'10\%', r'1\%', r'0.1\%']
colors = ['#81b56c', '#ec904e', '#de425b']

##8fb180
#e5dab2
#de9665

colors = np.flip(colors)
markers = ['*', 'x']

dumb_list = []

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(np.flip(false_alarm_prob_vec)):
 
	dumb = ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], color=colors[fa], label=r'$P_{\rm FA}=' + labels[fa] + '$')

	ax.plot(C_vec/L, hat_detected_ue_sig[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], marker=markers[0], markersize=6, color=colors[fa])
	ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], marker=markers[1], markersize=6, color=colors[fa])
	

	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	dumb_list.append(dumb)

dumb = ax.plot(C_vec/L, hat_detected_ue_sig[:, fa, :].mean(axis=-1), linewidth=0.0, linestyle='-', marker=markers[0], markersize=6, color='black', label='Signal-based')
dumb_list.append(dumb)
dumb = ax.plot(C_vec/L, hat_detected_ue_pow[:, fa, :].mean(axis=-1), linewidth=0.0, linestyle='--', marker=markers[1], markersize=6, color='black', label='Power-based')
dumb_list.append(dumb)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

ax.set_xlabel(r'Relative duration, $\frac{C}{L}$')
ax.set_ylabel(r'Probability of detection, $P_{\rm D}$')

ax.legend(fontsize='small', framealpha=0.5)

ax.set_xticks(xticks)

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5, alpha=0.5)

[(dumb.pop(0)).remove() for dumb in dumb_list]

#plt.tight_layout(pad=0.5)

plt.subplots_adjust(
	left = 0.18,  	# the left side of the subplots of the figure
	right = 0.99,   	# the right side of the subplots of the figure
	bottom = 0.15,  	# the bottom of the subplots of the figure
	top = 0.99,     	# the top of the subplots of the figure
	wspace = 0.5,  	# the amount of width reserved for space between subplots,
    	           	# expressed as a fraction of the average axis width
	hspace = 0.05   	# the amount of height reserved for space between subplots,
              	 	# expressed as a fraction of the average axis height
              )

plt.savefig('figures/detection_100x100.pdf', dpi='figure', format='pdf', transparent='True')


fig_size = (6.5/2, 3.5)
#xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
xticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]


#####


fig, ax = plt.subplots(figsize=fig_size)

dumb_list = []

# Go through all false alarm probabilities
for fa, false_alarm_prob in enumerate(np.flip(false_alarm_prob_vec)):

 
	#dumb = ax.plot(C_vec/L, distance_pow[:, fa, :].mean(axis=-1), linewidth=0.0, marker=markers[fa], markersize=6, color='black', label=r'$P_{\rm FA}=' + labels[fa] + '$')

	#avg_distance_sig = distance_sig[:, fa, :].mean(axis=-1)
	#avg_distance_sig[avg_distance_sig == 0.0] = np.nan

	#ax.plot(C_vec/L, avg_distance_sig, linewidth=1.5, linestyle=linestyles[fa], marker=markers[0], markersize=6, color=colors[fa])
	ax.plot(C_vec/L, distance_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle=linestyles[fa], marker=markers[1], markersize=6, color=colors[fa])
	

	#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, fa, :], axis=-1), np.max(mse_cha_est_pow[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
	#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, fa, :], axis=-1), np.max(mse_cha_est_sig[:, fa, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	#dumb_list.append(dumb)

# dumb = ax.plot(C_vec/L, distance_sig[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle='-', color='tab:blue', label='Signal-based')
# dumb_list.append(dumb)
# dumb = ax.plot(C_vec/L, distance_pow[:, fa, :].mean(axis=-1), linewidth=1.5, linestyle='--', color='tab:orange', label='Power-based')
# dumb_list.append(dumb)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

ax.set_xlabel(r'Relative duration, $\frac{C}{L}$')
ax.set_ylabel(r'Distance $\lVert \hat{\boldsymbol{\Theta}}_{\star} - \boldsymbol{\Theta}_{\star} \rVert^2_{{F}}$')

#ax.set_yscale('log')

#ax.legend(fontsize='x-small', framealpha=0.5)

#ax.set_ylim([10e-2, 0.25 * 10e-0])

ax.set_xticks(xticks)

ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5, alpha=0.5)

[(dumb.pop(0)).remove() for dumb in dumb_list]

plt.tight_layout(pad=0.5)

plt.subplots_adjust(
	left = 0.18,  	# the left side of the subplots of the figure
	right = 0.99,   	# the right side of the subplots of the figure
	bottom = 0.15,  	# the bottom of the subplots of the figure
	top = 0.99,     	# the top of the subplots of the figure
	wspace = 0.5,  	# the amount of width reserved for space between subplots,
    	           	# expressed as a fraction of the average axis width
	hspace = 0.05   	# the amount of height reserved for space between subplots,
              	 	# expressed as a fraction of the average axis height
              )

plt.savefig('figures/distance_100x100.pdf', dpi='figure', format='pdf', transparent='True')

#plt.show()


