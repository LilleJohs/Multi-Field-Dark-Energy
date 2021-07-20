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

w_list = np.linspace(start=-1, stop=-0.90, num=11)
omega_list = np.linspace(start=0.8, stop=1, num=11)
print(w_list)
print(omega_list)
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
        c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 1, gamma=4/3)
        c.run_background_eq_of_motion()
        #c.x_y_phase_plot()
        w_cur = c.get_eq_of_state()
        omega_cur = c.get_omega_phi()
        len_N = len(c.sol['t'])-1
        arrow_length = 1#sqrt((w_cur[len_N] - w)**2 + (omega_cur[len_N] - omega)**2)        
        
        arrow_x[j, i] = (w_cur[len_N] - w)/arrow_length
        arrow_y[j, i] = (omega_cur[len_N] - omega)/arrow_length

plt.quiver(w_list, omega_list, arrow_x, arrow_y)

plt.xlabel(r'$w_{\phi}$')
plt.ylabel(r'$\Omega_{\phi}$')
plt.title(r'$f(r) =r^{}$'.format(params['p']))

#plt.savefig('img/exp_spinning_p_-2_r_init.pdf', bbox_inches = 'tight')

plt.show()