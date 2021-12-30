import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.colors as mcolors

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')
plt.rc('font', size=15)
plt.rc('figure', autolayout=True)
plt.rc('axes', titlesize=16, labelsize=17)
plt.rc('lines', linewidth=2, markersize=6)
plt.rc('legend', fontsize=15)
plt.rc('text', usetex=True)

params = {
    'V0': 2.186,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 2*1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'potential': 'spinning',
    'metric': 'r_p'
}

x_range = [0.05, 1, 100]
colormap = ListedColormap(sns.cubehelix_palette(10, reverse = False).as_hex())
normalize = mcolors.Normalize(vmin=np.log(np.min(x_range)/10), vmax=np.log(np.max(x_range)*10))

colors = ['red', 'blue', 'green']
p_range=[0, 1, 2, 3]

fig, axs = plt.subplots(2, 2)
for m, p in enumerate(p_range):
    #print(np.floor(m/2), m%2)
    for i, x in enumerate(x_range):
        color = colormap(normalize(np.log(x_range[i])))
        params['r_init_multiplier'] = x
        params['p'] = p
        print(p, x)
        N_max = 12
        if p==2: N_max = 50
        elif p==1: N_max = 70
        c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = N_max, gamma=1)
        c.run_background_eq_of_motion()

        print('De Sitter Bound Lowest', min(c.get_de_sitter_bound()))
        field_derivative, delta_phi = c.get_field_derivative()
        print('Delta phi:', delta_phi)
        size = len(c.get_eq_of_state())
        w = c.get_eq_of_state()
        omega = c.get_omega_phi()
        N = c.sol['t']
        delta_phi_n = -1
        for j in range(size):
            cur = np.trapz(np.sqrt(3)*np.sqrt((1+w[:j])*omega[:j]), N[:j])
            if cur > 1:
                delta_phi_n = j
                break
        axs[int(np.floor(m/2)), m%2].plot(w, omega, label=r"$r_i = {{{}}}r_0$".format(x) +'\n dSB: '+"{:.2f}".format(min(c.get_de_sitter_bound())), color=color)
        if delta_phi_n>0: axs[int(np.floor(m/2)), m%2].plot(w[delta_phi_n], omega[delta_phi_n], 'go', color=color, linewidth=5, markersize=14)
 
    #cur_uni = Ellipse(xy=(-1, 0.7), width=0.065, height=0.02, 
                            #edgecolor='black', fc='black', lw=2)#, zorder=100)
 
    axs[int(np.floor(m/2)), m%2].set_xlim([-1.01, -0.8])
    axs[int(np.floor(m/2)), m%2].set_ylim([-0.02, 1.02])
    axs[int(np.floor(m/2)), m%2].set_title(r'$f(r) =r^{{{}}}$'.format(params['p']))
    #axs[int(np.floor(m/2)), m%2].axvline(x=-1, color='k', alpha=0.25)
   
    axs[int(np.floor(m/2)), m%2].legend()

axs[1, 0].set_xlabel(r'$w_{\phi}$')
axs[1, 1].set_xlabel(r'$w_{\phi}$')
axs[0, 0].set_ylabel(r'$\Omega_{\phi}$')
axs[1, 0].set_ylabel(r'$\Omega_{\phi}$')
fig.set_size_inches(7, 7)

plt.show()