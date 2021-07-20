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
    'r_init_multiplier': 1,
    'cosmo_constant': 0,
    'f0': 1
}

beta_ranges = [500, 1000, 1500]
for i, beta in enumerate(beta_ranges):
    params['beta'] = beta
    c = MultiFieldDarkEnergy(metric='exp', potential='exp_spinning', params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()
    N = c.sol['t']-8
    M_eff_s, cs_s = c.get_M_eff_squared()

    plt.plot(N, cs_s, label=r'$\beta = {{{}}}$'.format(beta))
plt.xlabel(r'$N$')
plt.ylabel(r'$c_s^2$')
plt.legend()
plt.grid()
plt.show()