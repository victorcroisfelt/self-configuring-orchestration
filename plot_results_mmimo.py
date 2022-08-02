import scipy.io

import numpy as np

import matplotlib.pyplot as plt

# LaTeX type definitions
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amssymb,amsmath,amsfonts,amsthm,mathtools,cuted,bbold} \usepackage[cmintegrals]{newtxmath}')

plt.rcParams['axes.titlesize'] = 12

########################################
# Load data
########################################

# Vector with sizes
sizes = ['50', '100', '250']

filenames = [
	'data/RESULTS_M64_N16_K16_L128_50x50.mat',
	'data/RESULTS_M64_N16_K16_L128_100x100.mat',
	'data/RESULTS_M64_N16_K16_L128_250x250.mat',
	#'data/RESULTS_M64_N16_K16_L128_500x500.mat'
]

# Dictionary with results type
type_res = {'mse': None, 'sinr': None, 'se': None}

# Create a dictionary to store the results for each enviroment
results = {
	'50': {'pow':type_res.copy(), 'sig':type_res.copy()},
	'100': {'pow':type_res.copy(), 'sig':type_res.copy()},
	'250': {'pow':type_res.copy(), 'sig':type_res.copy()},
 	#'500': {'pow': type_res.copy(), 'sig': type_res.copy()}
 	}

# Go through all sizes 
for ss, size in enumerate(sizes):

	mat = scipy.io.loadmat(filenames[ss])

	if ss == 0: 

		# Loading parameters 
		L = np.squeeze(mat['L'])
		C_vec = np.squeeze(mat['C_vec'])

	# Loading results 
	mse_cha_est_pow = mat['MSE_cha_est_pow']
	mse_cha_est_sig = mat['MSE_cha_est_sig']

	# sinr_bound_pow = mat['SINR_bound_pow']
	# sinr_bound_sig = mat['SINR_bound_sig']

	sinr_est_pow = mat['SINR_est_pow']
	sinr_est_sig = mat['SINR_est_sig']

	# se_bound_pow = mat['SE_bound_pow']
	# se_bound_sig = mat['SE_bound_sig']

	se_est_pow = mat['SE_est_pow']
	se_est_sig = mat['SE_est_sig']

	# Flatten results in terms of setup and channel realizations
	results[size]['pow']['mse'] = mse_cha_est_pow.reshape(*mse_cha_est_pow.shape[:2], -1)
	results[size]['sig']['mse'] = mse_cha_est_sig.reshape(*mse_cha_est_sig.shape[:2], -1)

	# sinr_bound_pow = sinr_bound_pow.reshape(*sinr_bound_pow.shape[:3], -1).mean(axis=0)
	# sinr_bound_sig = sinr_bound_sig.reshape(*sinr_bound_sig.shape[:3], -1).mean(axis=0)

	results[size]['pow']['sinr'] = sinr_est_pow.reshape(*sinr_est_pow.shape[:3], -1)
	results[size]['sig']['sinr'] = sinr_est_sig.reshape(*sinr_est_sig.shape[:3], -1)

	# se_bound_pow = se_bound_pow.reshape(*se_bound_pow.shape[:3], -1).mean(axis=0) #int(mat['tau_com'])/int(mat['tau_c']) * np.log2(sinr_est_pow) #
	# se_bound_sig = se_bound_sig.reshape(*se_bound_sig.shape[:3], -1).mean(axis=0) #int(mat['tau_com'])/int(mat['tau_c']) * np.log2(sinr_est_sig) #

	results[size]['pow']['se'] = se_est_pow.reshape(*se_est_pow.shape[:3], -1)
	results[size]['sig']['se'] = se_est_sig.reshape(*se_est_sig.shape[:3], -1)

########################################
# Plot
########################################
fig_size = (6.5, 3.5)
#xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
xticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
#xticks = [0, 0.3, 0.5, 0.7, 1.0]

linestyles = ['-', '--', '-.']
colors = ["#003f5c", "#bc5090", "#ffa600"]
markers = ['*', 'x']

fig, axes = plt.subplots(figsize=fig_size, ncols=3)

dumb_list = []

#dumb = ax.plot(C_vec/L, mse_cha_est_pow[:, -1, :].mean(axis=-1), linewidth=0.0, color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

# Go through all sizes 
for ss, size in enumerate(sizes):

	dumb = axes[0].plot(C_vec/L, results[size]['pow']['mse'][:, -1, :].mean(axis=-1), linewidth=1.0, linestyle=linestyles[ss], color=colors[ss], label=('$S = ' + str(size) + r'\times ' + str(size) + r'$'))
	
	axes[0].plot(C_vec/L, results[size]['sig']['mse'][:, -1, :].mean(axis=-1), linewidth=1.0, linestyle=linestyles[ss], color=colors[ss], marker=markers[0], markersize=5)
	axes[0].plot(C_vec/L, results[size]['pow']['mse'][:, -1, :].mean(axis=-1), linewidth=1.0, linestyle=linestyles[ss], color=colors[ss], marker=markers[1], markersize=5)

#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, -1, :], axis=-1), np.max(mse_cha_est_pow[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, -1, :], axis=-1), np.max(mse_cha_est_sig[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

	dumb_list.append(dumb[0])

dumb = axes[0].plot(C_vec/L, results[size]['sig']['mse'][:, -1, :].mean(axis=-1), linewidth=0.0, linestyle='-', marker=markers[0], markersize=5, color='black', label='Signal-based')
dumb_list.append(dumb[0])
dumb = axes[0].plot(C_vec/L, results[size]['sig']['mse'][:, -1, :].mean(axis=-1), linewidth=0.0, linestyle='--', marker=markers[1], markersize=5, color='black', label='Power-based')
dumb_list.append(dumb[0])

axes[0].set_title(r'Channel estimation')

axes[0].set_xlabel(r'Relative duration, $\frac{C}{L}$')
axes[0].set_ylabel('Average MSE per UE')

axes[0].set_yscale('log')

#axes[0].legend(fontsize='xx-small', framealpha=0.5)

axes[0].set_xticks(xticks)

axes[0].grid(color='#E9E9E9', linestyle='--', linewidth=0.5, alpha=0.5)



#plt.tight_layout()

#plt.savefig('figures/mse.pdf', dpi='figure', format='pdf', transparent='True')


#####


#fig, ax = plt.subplots(figsize=fig_size)


#dumb = ax.plot(C_vec/L, 10 * np.log10(sinr_est_pow[:, -1, :].mean(axis=-1)), linewidth=0.0, color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

# Go through all sizes 
for ss, size in enumerate(sizes):

	#ax.plot(C_vec/L, 10 * np.log10(np.median(results[size]['pow']['sinr'][:, :, -1, :], axis=0).mean(axis=-1)), linewidth=1.5, linestyle=linestyles[ss], marker=markers[0], markersize=5, color='tab:blue')
	#ax.plot(C_vec/L, 10 * np.log10(np.median(results[size]['sig']['sinr'][:, :, -1, :], axis=0).mean(axis=-1)), linewidth=1.5, linestyle=linestyles[ss], marker=markers[0], markersize=5, color='tab:orange')

	
	axes[1].plot(C_vec/L, 10 * np.log10(np.median(results[size]['sig']['sinr'][:, :, -1, :], axis=(-1, 0))), linewidth=1.5, linestyle=linestyles[ss], marker=markers[0], markersize=5, color=colors[ss])
	axes[1].plot(C_vec/L, 10 * np.log10(np.median(results[size]['pow']['sinr'][:, :, -1, :], axis=(-1, 0))), linewidth=1.5, linestyle=linestyles[ss], marker=markers[1], markersize=5, color=colors[ss])


#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, -1, :], axis=-1), np.max(mse_cha_est_pow[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, -1, :], axis=-1), np.max(mse_cha_est_sig[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

#ax.errorbar(C_vec/L, 10 * np.log10(sinr_est_sig[:, -1, :].mean(axis=-1)), 10 * np.log10(sinr_est_sig[:, -1, :].std(axis=-1)))

#dumb_list.append(dumb)

axes[1].set_title(r'Comm.: SINR')

axes[1].set_xlabel(r'Relative duration, $\frac{C}{L}$')
axes[1].set_ylabel('Median SINR per UE [dB]')

#axes[0].legend(fontsize='x-small', framealpha=0.5)

axes[1].set_xticks(xticks)

axes[1].grid(color='#E9E9E9', linestyle='--', linewidth=0.5, alpha=0.5)

#[(dumb.pop(0)).remove() for dumb in dumb_list]

# plt.tight_layout()

# plt.savefig('figures/sinr.pdf', dpi='figure', format='pdf', transparent='True')


######


# fig, ax = plt.subplots(figsize=fig_size)



#dumb = ax.plot(C_vec/L, 10 * np.log10(sinr_est_pow[:, -1, :].mean(axis=-1)), linewidth=0.0, color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

# Go through all sizes 
for ss, size in enumerate(sizes):

	axes[2].plot(C_vec/L, np.median(results[size]['sig']['se'][:, :, -1, :], axis=(-1, 0)), linewidth=1.5, linestyle=linestyles[ss], marker=markers[0], markersize=5, color=colors[ss])
	axes[2].plot(C_vec/L, np.median(results[size]['pow']['se'][:, :, -1, :], axis=(-1, 0)), linewidth=1.5, linestyle=linestyles[ss], marker=markers[1], markersize=5, color=colors[ss])

#ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, -1, :], axis=-1), np.max(mse_cha_est_pow[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
#ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, -1, :], axis=-1), np.max(mse_cha_est_sig[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

#ax.errorbar(C_vec/L, 10 * np.log10(sinr_est_sig[:, -1, :].mean(axis=-1)), 10 * np.log10(sinr_est_sig[:, -1, :].std(axis=-1)))

#dumb_list.append(dumb)

axes[2].set_title(r'Comm.: SE')

axes[2].set_xlabel(r'Relative duration, $\frac{C}{L}$')
axes[2].set_ylabel('Median SE per UE [bits/s/Hz]')

#axes[2].legend(fontsize='x-small', framealpha=0.5)

axes[2].set_xticks(xticks)

axes[2].grid(color='#E9E9E9', linestyle='--', linewidth=0.5, alpha=0.5)

#[(dumb.pop(0)).remove() for dumb in dumb_list]

for a in axes.flatten():
    a.tick_params(axis='both', which='major', labelsize=10)
    a.tick_params(axis='both', which='minor', labelsize=10)


#fig.subplots_adjust(bottom=0.3, wspace=0.33)

axes[1].legend(
			handles=dumb_list,
			fontsize='x-small', framealpha=0.5,
			loc='upper center', 
            bbox_to_anchor=(0.5, -0.3), fancybox=False, shadow=False, ncol=5
            )

[dumb.remove() for dumb in dumb_list]

plt.subplots_adjust(
	left = 0.1,  	# the left side of the subplots of the figure
	right = 0.97,   	# the right side of the subplots of the figure
	bottom = 0.28,  	# the bottom of the subplots of the figure
	top = 0.93,     	# the top of the subplots of the figure
	wspace = 0.5,  	# the amount of width reserved for space between subplots,
    	           	# expressed as a fraction of the average axis width
	hspace = 0.05   	# the amount of height reserved for space between subplots,
              	 	# expressed as a fraction of the average axis height
              )

plt.show()

#plt.savefig('figures/se.pdf', dpi='figure', format='pdf', transparent='True')


# #####


# fig, ax = plt.subplots(figsize=fig_size)

# dumb_list = []


# #dumb = ax.plot(C_vec/L, mse_cha_est_pow[:, -1, :].mean(axis=-1), linewidth=0.0, color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

# ax.plot(C_vec/L, se_est_sig[:, -1, :].mean(axis=-1), linewidth=1.5, linestyle='-', color='tab:blue')
# ax.plot(C_vec/L, se_est_pow[:, -1, :].mean(axis=-1), linewidth=1.5, linestyle='--', color='tab:orange')

# #ax.errorbar(C_vec/L, se_est_sig[:, -1, :].mean(axis=-1), np.std(se_est_sig[:, -1, :]))


# #ax.fill_between(C_vec/L, np.percentile(se_est_sig[:, -1, :], 25, axis=-1), np.percentile(se_est_sig[:, -1, :], 75, axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
# #ax.fill_between(C_vec/L, np.percentile(se_est_pow[:, -1, :], 25, axis=-1), np.percentile(se_est_pow[:, -1, :], 75, axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

# #dumb_list.append(dumb)

# ax.set_xlabel(r'Relative duration, $\frac{C}{L}$')
# ax.set_ylabel('SE [bit/s/Hz]')

# ax.legend(fontsize='x-small', framealpha=0.5)

# ax.set_xticks(xticks)

# ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

# [(dumb.pop(0)).remove() for dumb in dumb_list]

# plt.tight_layout()

#plt.show()


#####


# fig, ax = plt.subplots(figsize=fig_size)

# dumb_list = []


# #dumb = ax.plot(C_vec/L, 10 * np.log10(sinr_est_pow[:, -1, :].mean(axis=-1)), linewidth=0.0, color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

# ax.plot(C_vec/L, 10 * np.log10(sinr_bound_sig[:, -1, :].mean(axis=-1)), linewidth=1.5, linestyle='-', color='tab:blue')
# ax.plot(C_vec/L, 10 * np.log10(sinr_bound_pow[:, -1, :].mean(axis=-1)), linewidth=1.5, linestyle='--', color='tab:orange')

# #ax.fill_between(C_vec/L, np.min(mse_cha_est_pow[:, -1, :], axis=-1), np.max(mse_cha_est_pow[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
# #ax.fill_between(C_vec/L, np.min(mse_cha_est_sig[:, -1, :], axis=-1), np.max(mse_cha_est_sig[:, -1, :], axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

# #ax.errorbar(C_vec/L, 10 * np.log10(sinr_est_sig[:, -1, :].mean(axis=-1)), 10 * np.log10(sinr_est_sig[:, -1, :].std(axis=-1)))

# #dumb_list.append(dumb)

# ax.set_xlabel(r'Relative duration, $\frac{C}{L}$')
# ax.set_ylabel('SINR [dB]')

# #ax.legend(fontsize='x-small', framealpha=0.5)

# ax.set_xticks(xticks)

# ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

# [(dumb.pop(0)).remove() for dumb in dumb_list]

# plt.tight_layout()

# plt.savefig('figures/sinr.pdf', dpi='figure', format='pdf', transparent='True')


# #####


# fig, ax = plt.subplots(figsize=fig_size)

# dumb_list = []


# #dumb = ax.plot(C_vec/L, mse_cha_est_pow[:, -1, :].mean(axis=-1), linewidth=0.0, color='black', label=r'$P_{\rm FA}=' + str(false_alarm_prob) + '$')

# ax.plot(C_vec/L, se_bound_sig[:, -1, :].mean(axis=-1), linewidth=1.5, linestyle='-', color='tab:blue')
# ax.plot(C_vec/L, se_bound_pow[:, -1, :].mean(axis=-1), linewidth=1.5, linestyle='--', color='tab:orange')

# #ax.errorbar(C_vec/L, se_est_sig[:, -1, :].mean(axis=-1), np.std(se_est_sig[:, -1, :]))


# #ax.fill_between(C_vec/L, np.percentile(se_est_sig[:, -1, :], 25, axis=-1), np.percentile(se_est_sig[:, -1, :], 75, axis=-1), linewidth=0.0, alpha=0.5, color='tab:blue')
# #ax.fill_between(C_vec/L, np.percentile(se_est_pow[:, -1, :], 25, axis=-1), np.percentile(se_est_pow[:, -1, :], 75, axis=-1), linewidth=0.0, alpha=0.5, color='tab:orange')

# #dumb_list.append(dumb)

# ax.set_xlabel(r'Relative duration, $\frac{C}{L}$')
# ax.set_ylabel('SE [bit/s/Hz]')

# ax.legend(fontsize='x-small', framealpha=0.5)

# ax.set_xticks(xticks)

# ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

# [(dumb.pop(0)).remove() for dumb in dumb_list]

# plt.tight_layout()

# #plt.show()



