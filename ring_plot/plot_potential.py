import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
import matplotlib.pylab as plt
import numpy as np
from numpy import pi, cos, sin

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')
plt.rc('font', size=15)
plt.rc('figure', autolayout=True)
plt.rc('axes', titlesize=16, labelsize=17)
plt.rc('lines', linewidth=2, markersize=6)
plt.rc('legend', fontsize=15)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

params = {
    'V0': 2.15,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 5*1e-3,
    'x_p_init': 0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 4,
    'p': 2,
    'cosmo_constant': 0
}

c = MultiFieldDarkEnergy(metric='r_p', potential='exp_spinning', params=params, N_min = 0, N_max = 7.02, gamma=1)
r = np.linspace(0, 0.005, 100)
theta = np.linspace(0, 4*pi, 100)
R, THETA = np.meshgrid(r, theta)
Z, _, _ = c.get_V_and_diff(R, THETA)
ax.plot_surface(R*cos(THETA), R*sin(THETA), Z, alpha=0.3)

c.run_background_eq_of_motion()
phi = c.sol['y'][3]
theta = c.sol['y'][4]
#theta = theta * 2*pi / max(theta)
pot, _, _ = c.get_V_and_diff(phi, theta)
ax.plot(phi*cos(theta), phi*sin(theta), pot, label='parametric curve', linewidth=5)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zlabel(r'$V \; [H_0^2M_{Pl}^2]$')

plt.show()