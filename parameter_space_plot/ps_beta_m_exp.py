import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

params = {
    'V0': 2.15,
    'm': 80,
    'r0': 0,
    'alpha': np.sqrt(2),
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'beta': 100,
    'r_init': 0,
    'cosmo_constant': 0,
    'f0': 1
}

cur_time = Polygon([(-1.05, 0.65), (-0.95, 0.75), (-0.95, 0.65), (-1.05, 0.75)])

def get_sol(params):
    c = MultiFieldDarkEnergy(metric='exp', potential='exp_spinning', params=params, N_min = 0, N_max = 9, gamma=1)
    c.run_background_eq_of_motion()
    #c.plot_swampland_bound()
    #field_derivative, delta_phi = c.get_field_derivative()
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
                #cur = np.trapz(np.sqrt(3)*np.sqrt((1+w[:k])*omega[:k]), N[:k])
                # cur = delta phi at omega=0.7 w=-1
                cur_param = 0.5
                #if cur > 1:
                #    cur_param= 1
                break
            elif k == size-1:
                # This solution has NEVER been in omega=0.7 w=-1
                cur_param = 1
    if omega[size-1] < 0.9 and w[size-1] < -0.8:
        cur_param=-1
    

    return cur_param, np.nanmin(c.get_M_eff_squared())

length = 50
list_accepted = np.zeros((length, length))

beta_range = np.linspace(0, 2000, length)
m_range = np.linspace(0, 300, length)
colors = ['red', 'blue']

meff_over_beta = np.zeros(length)

for i, beta in enumerate(beta_range):
    print(i)
    for j, m in enumerate(m_range):
        params['beta'] = beta
        params['m'] = m
        
        cur_param = get_sol(params)

        list_accepted[i, j], meff = cur_param
        print(meff)
        if meff > 0 and meff_over_beta[i] == 0:
           
            meff_over_beta[i] = m
        

fig, axs = plt.subplots(1)
axs.set_xlabel(r'$m [H_0]$')
axs.set_ylabel(r'$\beta$')
#axs.plot(meff_over_beta[:np.argmin(meff_over_beta)], beta_range[:np.argmin(meff_over_beta)])
plt.fill_between(meff_over_beta[:np.argmin(meff_over_beta)], beta_range[:np.argmin(meff_over_beta)])
axs.imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(beta_range), np.amax(beta_range)], aspect='auto',cmap='RdGy')

plt.show()