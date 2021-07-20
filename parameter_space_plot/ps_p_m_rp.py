import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

params = {
    'Omega_M_0': 0.3,
    'Omega_R_0': 6e-5,
    'V0': 2.15,
    'm': 10,
    'r0': 7*1e-4,
    'alpha': 1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'p': 2,
    'r_init_multiplier': 1,
    'cosmo_constant': 0,
}

cur_time = Polygon([(-1.05, 0.65), (-0.95, 0.75), (-0.95, 0.65), (-1.05, 0.75)])

def get_sol(params):
    c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 9, gamma=1)
    c.run_background_eq_of_motion()
    #c.plot_swampland_bound()
    field_derivative, delta_phi = c.get_field_derivative()
    size = len(c.get_eq_of_state())
    w = c.get_eq_of_state()
    omega = c.get_omega_phi()
    N = c.sol['t']
    cur_param = 0.5
    if min(c.get_de_sitter_bound()) < 0.5:
        cur_param = 0
        #print(p, m)
        #print('De Sitter violated', min(c.get_de_sitter_bound()))
    else:
        for k in range(size):
            point = Point(w[k], omega[k])
            if cur_time.contains(point):
                # This solution has once been in omega=0.7 w=-1
                #print('H:', c.get_H()[k])
                cur = np.trapz(np.sqrt(3)*np.sqrt((1+w[:k])*omega[:k]), N[:k])
                # cur = delta phi at omega=0.7 w=-1
                cur_param = 0.5 + cur
                if cur > 1:
                    cur_param= 1
                break
            elif k == size-1:
                # This solution has NEVER been in omega=0.7 w=-1
                cur_param = 1
    if omega[size-1] < 0.9 and w[size-1] < -0.8:
        cur_param=-1  
    return cur_param

length = 10
list_accepted = np.zeros((length, length))

p_range = np.linspace(1.6, 3, length)
m_range = np.linspace(0, 400, length)
colors = ['red', 'blue']

for i, p in enumerate(p_range):
    print(i)
    for j, m in enumerate(m_range):
        params['p'] = p
        params['m'] = m
        
        cur_param = get_sol(params)

        list_accepted[i, j] = cur_param
        #if cur_param == 1:
        #    print('Last eq of state:', w[size-1])
        #print(cur_param, p, m)

fig, axs = plt.subplots(2)
axs[0].set_xlabel(r'$m [H_0]$')
axs[0].set_ylabel(r'$p$')
axs[0].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(p_range), np.amax(p_range)], aspect='auto',cmap='RdGy')

params['p'] = 2
alpha_range = np.linspace(0, 0.003, length)
m_range = np.linspace(0, 400, length)
colors = ['red', 'blue']

for i, alpha in enumerate(alpha_range):
    print(i)
    for j, m in enumerate(m_range):
        params['alpha'] = alpha
        params['m'] = m
        
        cur_param = get_sol(params)

        list_accepted[i, j] = cur_param
        #if cur_param == 1:
        #    print('Last eq of state:', w[size-1])
        #print(cur_param, p, m)

axs[1].set_xlabel(r'$m [H_0]$')
axs[1].set_ylabel(r'$\alpha$')
axs[1].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto',cmap='RdGy')

#plt.savefig('img/test.pdf', bbox_inches = 'tight')

plt.show()