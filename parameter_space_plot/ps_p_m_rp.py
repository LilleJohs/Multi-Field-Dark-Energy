import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')
plt.rc('font', size=15)
plt.rc('figure', autolayout=True)
plt.rc('axes', titlesize=16, labelsize=17)
plt.rc('lines', linewidth=2, markersize=6)
plt.rc('legend', fontsize=15)

params = {
    'V0': 2.15,
    'm': 50,
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
                #cur = np.trapz(np.sqrt(3)*np.sqrt((1+w[:k])*omega[:k]), N[:k])
                # cur = delta phi at omega=0.7 w=-1
                cur_param = 0.5# + cur
                break
            elif k == size-1:
                # This solution has NEVER been in omega=0.7 w=-1
                cur_param = 1
    if omega[size-1] < 0.9 and w[size-1] < -0.8:
        cur_param=-1  
    return cur_param, np.nanmin(c.get_M_eff_squared())

length = 40
list_accepted = np.zeros((length, length))

p_range = np.linspace(1.6, 3.5, length)
m_range = np.linspace(0, 500, length)
meff_over_p = np.zeros(length)
colors = ['red', 'blue']

for j, m in enumerate(m_range):
    print(j)
    for i, p in enumerate(p_range):
        params['p'] = p
        params['m'] = m
        
        cur_param = get_sol(params)

        list_accepted[i, j], meff = cur_param
        #print(meff)
        if meff < 0 and meff_over_p[j] == 0:
            meff_over_p[j] = p
        elif i == length-1 and meff_over_p[j] == 0:
            meff_over_p[j] = p

fig, axs = plt.subplots(2)
axs[0].set_xlabel(r'$m [H_0]$')
axs[0].set_ylabel(r'$p$')
axs[0].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(p_range), np.amax(p_range)], aspect='auto',cmap='RdGy')
axs[0].fill_between(m_range, meff_over_p, y2 = p_range[length-1], color='blue', alpha=0.4)

params['p'] = 2
alpha_range = np.linspace(0, 0.015, length)
m_range = np.linspace(0, 500, length)
meff_over_alpha = np.zeros(length)
colors = ['red', 'blue']

for j, m in enumerate(m_range):
    print(j)
    for i, alpha in enumerate(alpha_range):
        params['alpha'] = alpha
        params['m'] = m
        
        cur_param = get_sol(params)

        list_accepted[i, j], meff = cur_param
        #print(meff)
        if meff < 0 and meff_over_alpha[j] == 0:
            meff_over_alpha[j] = alpha
        elif i == length-1 and meff_over_alpha[j] == 0:
            meff_over_alpha[j] = alpha
axs[1].set_xlabel(r'$m [H_0]$')
axs[1].set_ylabel(r'$\alpha$')
axs[1].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto',cmap='RdGy')
axs[1].fill_between(m_range, meff_over_alpha, y2 = alpha_range[length-1], color='blue', alpha=0.4)

#plt.savefig('img/test.pdf', bbox_inches = 'tight')
fig.set_size_inches(7, 7)

plt.show()