import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

params = {
    'V0': 2.15,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 5*1e-3,
    'x_p_init': 0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'cosmo_constant': 0
}

r_init_ranges = [0.2, 8]
for i, r_init in enumerate(r_init_ranges):
    params['r_init_multiplier'] = r_init
    c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()
    phi = c.sol['y'][3]
    theta = c.sol['y'][4]

    p = params['p']
    req = np.power(p*params['alpha']**2 * 2.15*np.exp(-params['alpha']*theta)/(6*params['m']**2), 1/(2+p))

    plt.polar(theta/30, phi)
plt.polar(np.linspace(0, 2*pi, 100), np.ones(100)*params['r0'], label=r'$r_0$')
#plt.polar(np.linspace(0, 2*pi, 100), req*np.ones(100), label=r'$r_{eq}$', c='k', linestyle='dashed', linewidth=3)
plt.polar(theta/30, req, label=r'$r_{eq}$', c='k', linestyle='dashed', linewidth=3)
plt.legend()
plt.grid(True)
plt.box(on=None)
plt.xticks([])
plt.show()