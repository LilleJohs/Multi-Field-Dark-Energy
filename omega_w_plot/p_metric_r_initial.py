import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from matplotlib.patches import Ellipse

params = {
    'Omega_M_0': 0.3,
    'Omega_R_0': 6e-5,
    'V0': 2.186,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'cosmo_constant': 0
}

x_range = [0.01, 1, 100]
colors = ['red', 'blue', 'green']
p_range=[0, 1, 2, 3]

fig, axs = plt.subplots(2, 2)
for m, p in enumerate(p_range):
    #print(np.floor(m/2), m%2)
    for i, x in enumerate(x_range):
        params['r_init_multiplier'] = x
        params['p'] = p
        print(p, x)
        N_max = 12
        if p==2: N_max=50
        c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = N_max, gamma=1)
        c.run_background_eq_of_motion()
        #c.x_y_phase_plot()
        #c.plot_swampland_bound()

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
        axs[int(np.floor(m/2)), m%2].plot(w, omega, label=r"$r_i = {{{}}}r_0$".format(x) +'\n dSB: '+"{:.2f}".format(min(c.get_de_sitter_bound())), color=colors[i])
        if delta_phi_n>0: axs[int(np.floor(m/2)), m%2].plot(w[delta_phi_n], omega[delta_phi_n], 'go', color=colors[i], linewidth=5, markersize=14)
    #plt.xlim([-1.02,  1.02])
    #plt.ylim([-0.05, 1.05])
    #cur_uni = Ellipse(xy=(-1, 0.7), width=0.065, height=0.02, 
                            #edgecolor='black', fc='black', lw=2)#, zorder=100)
    #ax = plt.gca()
    #ax.add_patch(cur_uni)
    #axs[int(np.floor(m/2)), m%2].xlabel(r'$w_{\phi}$')
    #axs[int(np.floor(m/2)), m%2].ylabel(r'$\Omega_{\phi}$')
    axs[int(np.floor(m/2)), m%2].set_xlim([-1.01, -0.8])
    axs[int(np.floor(m/2)), m%2].set_ylim([-0.02, 1.02])
    axs[int(np.floor(m/2)), m%2].set_title(r'$f(r) =r^{{{}}}$'.format(params['p']))
    
   
    axs[int(np.floor(m/2)), m%2].legend()

axs[1, 0].set_xlabel(r'$w_{\phi}$')
axs[1, 1].set_xlabel(r'$w_{\phi}$')
axs[0, 0].set_ylabel(r'$\Omega_{\phi}$')
axs[1, 0].set_ylabel(r'$\Omega_{\phi}$')
#plt.savefig('p_init.pdf', bbox_inches = 'tight')
plt.show()