from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from matplotlib.patches import Ellipse

params = {
    'Omega_M_0': 0.3,
    'Omega_R_0': 6e-5,
    'V0': 2.186,
    'm': 80,
    'r0': 7*1e-4,
    'alpha': 9.2e-4,
    'x_p_init': 0.001,
    'x_t_init': 0.001,
    'y_1_init': 0.001,
    'r_init_multiplier': 1,
    'p': 3,
    'cosmo_constant': 0
}

x_range = [0.0007]
colors = ['red', 'blue', 'green']

for i, x in enumerate(x_range):
    #omega=0.9
    #w=-0.99
    #params['y_1_init'] = sqrt(omega*(1-w)/2)
    #params['x_t_init'] = sqrt(omega*(1+w)/2)#/sqrt(2)
    #params['x_p_init'] = sqrt(omega*(1+w)/2)#/sqrt(2)

    #params['V0'] = V0
    c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 3, gamma=1)
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
    print(w[delta_phi_n], omega[delta_phi_n])
    plt.plot(w, omega, label=r"$x_r = x_{\theta} = $"+str(x) +'\n dSB: '+"{:.2f}".format(min(c.get_de_sitter_bound())), color=colors[i])
    if delta_phi_n>0: plt.plot(w[delta_phi_n], omega[delta_phi_n], 'go', color='black')
#plt.xlim([-1.02,  1.02])
#plt.ylim([-0.05, 1.05])
#cur_uni = Ellipse(xy=(-1, 0.7), width=0.065, height=0.02, 
                        #edgecolor='black', fc='black', lw=2)#, zorder=100)
#ax = plt.gca()
#ax.add_patch(cur_uni)
plt.xlabel(r'$w_{\phi}$')
plt.ylabel(r'$\Omega_{\phi}$')
plt.title(r'$f(r) =r^{}$'.format(params['p']))
plt.legend()
#plt.savefig('img/exp_spinning_p_-2_r_init.pdf', bbox_inches = 'tight')

plt.show()