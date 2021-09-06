import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from scipy.optimize import fsolve

params = {
    'V0': 2.15,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 1e-3,
    'x_p_init': 0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'cosmo_constant': 0
}
params = {
    'V0': 2.186,
    'm': 50,
    'r0': 0,
    'alpha': np.sqrt(2),
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'beta': 1000,
    'r_init': 0,
    'cosmo_constant': 0,
    'f0': 1
}

def solve_for_req(params, theta):
    func = lambda req : req**3 + 2*params['V0']/params['m']**2 * np.exp(-params['alpha'] * theta) * req - (params['V0']*params['alpha']*np.exp(-params['alpha']*theta)/(params['m']**2))**2 * params['beta']/3 * np.exp(-params['beta']*req)
    sol = fsolve(func, 0.01)
    return sol

r_init_ranges = [0]
for i, r_init in enumerate(r_init_ranges):
    params['r_init'] = r_init
    c = MultiFieldDarkEnergy(metric='exp', potential='exp_spinning', params=params, N_min = 0, N_max = 20, gamma=1)
    c.run_background_eq_of_motion()
    phi = c.sol['y'][3]
    theta = c.sol['y'][4]
    #p=params['p']
    #req = np.power(p*params['alpha']**2 * 2.15*np.exp(-params['alpha']*theta)/(6*params['m']**2), 1/(2+p))
    req = np.zeros(len(theta))
    for i in range(len(theta)):
        req[i] = solve_for_req(params, theta[i])
    plt.polar(10*theta, phi)
#plt.polar(np.linspace(0, 2*pi, 100), np.ones(100)*params['r0'], label=r'$r_0$')
#plt.polar(np.linspace(0, 2*pi, 100), req*np.ones(100), label=r'$r_{eq}$', c='k', linestyle='dashed', linewidth=3)
plt.polar(10*theta, req, label=r'$r_{eq}$', c='k', linestyle='dashed', linewidth=3)
plt.legend()
plt.grid(True)
plt.box(on=None)
plt.xticks([])
plt.show()