import scipy.io



import numpy as np



import matplotlib.pyplot as plt

from matplotlib import rc



# LaTeX type definitions

rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern'], 'size': 14})

#rc('text', usetex=True)



#plt.rcParams.update({

#    "text.usetex": True,

#    "font.family": "sans-serif",

#    "font.sans-serif": ["Helvetica"]})





########################################

# Load data

########################################

mat = scipy.io.loadmat('RES_Figure_effect_M64_N16_K16_L32_250x250.mat')



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

ax.set_ylabel(r'$\|\|\mathring{\mathbf{h}}_k\|\|_2^2$ [dB]')





ax.set_xticks([8, 16, 24, 32, 40, 48, 64])

#ax.grid(color='#E9E9E9', linestyle='--', linewidth=0.5)



plt.tight_layout()

plt.savefig('figures/hris_effect.pdf', dpi='figure', format='pdf', transparent='True')



#####





