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
    'm': 80,
    'r0': 7*1e-4,
    'alpha': 1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'cosmo_constant': 0
}

w_eq = -1 + 1/(1+sqrt(3*params['V0'])/(params['m']*params['alpha']))
print(w_eq)

def make_arrows(w_list, omega_list):
    arrow_x = np.zeros((len(w_list), len(omega_list)))
    arrow_y = np.zeros((len(w_list), len(omega_list)))

    for i, w in enumerate(w_list):
        print(i)
        for j, omega in enumerate(omega_list):
            # y^2 = omega*(1-w)/2
            # x^2 = omega*(1+w)/2
            params['y_1_init'] = sqrt(omega*(1-w)/2)
            #params['x_t_init'] = sqrt(omega*(1+w)/2)#/sqrt(2)
            params['x_p_init'] = sqrt(omega*(1+w)/2)#/sqrt(2)
            c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 2, gamma=4/3)

            c.run_background_eq_of_motion()
            #c.x_y_phase_plot()
            w_cur = c.get_eq_of_state()
            omega_cur = c.get_omega_phi()
            len_N = len(c.sol['t'])-1
            arrow_length = 1#sqrt((w_cur[len_N] - w)**2 + (omega_cur[len_N] - omega)**2)        
            
            arrow_x[j, i] = (w_cur[len_N] - w)/arrow_length
            arrow_y[j, i] = (omega_cur[len_N] - omega)/arrow_length
    return arrow_x, arrow_y

num = 11
w_list = np.linspace(start=-0.999, stop=-0.9, num=num)
omega_list = np.linspace(start=0.9, stop=0.9999, num=num)
arrow_x, arrow_y = make_arrows(w_list, omega_list)

red_box_x = [-1, -0.9, -0.9, -1, -1]
red_box_y = [0.9, 0.9, 1, 1, 0.9]

fig, axs = plt.subplots(2)
axs[0].set_xlabel(r'$w_{\phi}$')
axs[0].set_ylabel(r'$\Omega_{\phi}$')
axs[0].quiver(w_list, omega_list, arrow_x, arrow_y)
axs[0].plot(w_eq, 1, 'go', color='blue')

w_list = np.linspace(start=-0.99, stop=0.95, num=num)
omega_list = np.linspace(start=0.01, stop=0.99, num=num)
arrow_x, arrow_y = make_arrows(w_list, omega_list)
axs[1].set_xlabel(r'$w_{\phi}$')
axs[1].set_ylabel(r'$\Omega_{\phi}$')
axs[1].quiver(w_list, omega_list, arrow_x, arrow_y)
axs[1].plot(red_box_x, red_box_y, c='r')
#axs[1].plot(w_eq, 1, 'go', color='blue')

#plt.savefig('p_init_arrow.pdf', bbox_inches = 'tight')

plt.show()