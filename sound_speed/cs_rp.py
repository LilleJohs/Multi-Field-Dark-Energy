import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

params = {
    'V0': 2.15,
    'm': 100,
    'r0': 7*1e-4,
    'alpha': 1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'cosmo_constant': 0
}

p_ranges = [0, 1, 2]
for i, p in enumerate(p_ranges):
    params['p'] = p
    c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()
    N = c.sol['t']-8
    Vnn, omega = c.get_turning_rate()
    M_eff_s, cs_s = c.get_M_eff_squared()
    '''
    plt.plot(N, omega**2, label=r'$\Omega^2$')
    plt.plot(N, M_eff_s, label=r'$M_{eff}^2$')
    plt.legend()
    plt.figure()'''
    #plt.plot(N, cs_s, label=r'$c_s^2$')
    #plt.plot(N, 1+4*omega**2/M_eff_s, label=r'$1+4\Omega^2/M_{eff}^2$')
    plt.plot(N, cs_s, label=r'$p = {{{}}}$'.format(p))
plt.xlabel(r'$N$')
plt.ylabel(r'$c_s^2$')
plt.legend()
plt.grid()
plt.show()