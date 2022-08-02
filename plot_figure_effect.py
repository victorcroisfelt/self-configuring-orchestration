import scipy.io



import numpy as np



import matplotlib.pyplot as plt

from matplotlib import rc



# LaTeX type definitions
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amssymb,amsmath,amsfonts,amsthm,mathtools,cuted,bbold} \usepackage[cmintegrals]{newtxmath}')

#plt.rcParams.update({

#    "text.usetex": True,

#    "font.family": "sans-serif",

#    "font.sans-serif": ["Helvetica"]})





########################################

# Load data

########################################

mat = scipy.io.loadmat('data/RES_Figure_effect_M64_N16_K16_L32_250x250.mat')



# Loading parameters 

L = np.squeeze(mat['L'])

C_vec = np.squeeze(mat['C_vec'])

false_alarm_prob_vec = np.squeeze(mat['false_alarm_prob_vec'])



# Loading results 

channel_power_no___hris = mat['channel_power_no_hris']

channel_power_with_hris = mat['channel_power_with_hris']







########################################

# Plot

########################################

fig_size = (6.5, 2.5)



fig, ax = plt.subplots(figsize=fig_size)

colors = [[248/255, 203/255, 173/255], [255/255, 242/255, 204/255]]





markers = ['.', '+', 'x']

dumb_list = []



duration = np.array(range(channel_power_with_hris.shape[1]))



# Go through all false alarm probabilities

ax.plot(duration+1, 10*np.log10(channel_power_with_hris[0,:]), linewidth=1.5, linestyle='-',  color='red' , label='w/  HRIS')

ax.plot(duration+1, 10*np.log10(channel_power_no___hris[0,:]), linewidth=1.5, linestyle='--', color='black', label='w/o HRIS')





style = dict(size=14, color='black', alpha=1, ha='center', va='center')

ax.text(16, -52, 'Probing', fontdict=None, **style)

ax.text(48, -52, 'Reflection', fontdict=None, **style)



ax.text(4, -64, r'$\mathring{\mathbf{h}}_{kt}$', fontdict=None, **style)

ax.text(34, -64, r'$\mathring{\mathbf{h}}_{k}^{\star}$', fontdict=None, **style)





ax.set_xlim([1, 63])



#ax.arrow(x=0, dx=10, y=-50, dy=0, color='gray')



#ax.axvspan(1, 32, facecolor='orange', alpha=0.5)

#ax.axvspan(32, 64, facecolor='yellow', alpha=0.5)

ax.axvspan(1, 32, facecolor=colors[0], alpha=1)

ax.axvspan(32, 64, facecolor=colors[1], alpha=1)



ax.legend()



ax.set_xlabel(r'Subblock $t$')

ax.set_ylabel('Intensity of equivalent \n channel, $\lVert\mathring{\mathbf{h}}_k\lVert_2^2$ [dB]')





ax.set_xticks([1, 8, 16, 24, 32, 40, 48, 56, 64])

#ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.subplots_adjust(
	left = 0.14,  	# the left side of the subplots of the figure
	right = 0.99,   	# the right side of the subplots of the figure
	bottom = 0.17,  	# the bottom of the subplots of the figure
	top = 0.82,     	# the top of the subplots of the figure
	wspace = 0.5,  	# the amount of width reserved for space between subplots,
    	           	# expressed as a fraction of the average axis width
	hspace = 0.05   	# the amount of height reserved for space between subplots,
              	 	# expressed as a fraction of the average axis height
              )


plt.savefig('figures/hris_effect.pdf', dpi='figure', format='pdf', transparent='True')

plt.show()

#####





