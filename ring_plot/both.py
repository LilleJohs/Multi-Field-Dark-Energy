import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
import matplotlib.pylab as plt
import numpy as np
from numpy import pi, cos, sin
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm

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
ax = fig.add_subplot(1, 2, 2, projection='3d')
r_init_ranges = [0.2, 8]
colormap = ListedColormap(sns.cubehelix_palette(16, reverse = False).as_hex())
normalize = mcolors.Normalize(vmin=min(r_init_ranges), vmax=max(r_init_ranges)+15)

params = {
    'V0': 2.186,
    'm': 200,
    'r0': 7*1e-4,
    'alpha': 2*1e-3,
    'x_p_init': 0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 4,
    'p': 2,
    'metric': 'r_p',
    'potential': 'spinning'
}

c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 7.17, gamma=1)
r = np.linspace(0, 0.003, 100)
theta = np.linspace(0, 4*pi, 100)
R, THETA = np.meshgrid(r, theta)
Z, _, _ = c.get_V_and_diff(R, THETA)
ax.plot_surface(R*cos(THETA), R*sin(THETA), Z, alpha=0.5, color=colormap(normalize(10*r_init_ranges[0])))

c.run_background_eq_of_motion()
phi = c.sol['y'][3]
theta = c.sol['y'][4]
#theta = theta * 2*pi / max(theta)
pot, _, _ = c.get_V_and_diff(phi, theta)
ax.plot(phi*cos(theta), phi*sin(theta), pot, label='parametric curve', linewidth=5, color=colormap(normalize(r_init_ranges[1])))
ax.set_xticklabels([])
ax.set_yticklabels([])


ax.set_zlabel(r'$V \; [H_0^2]$', labelpad=10)

ax = fig.add_subplot(1, 2, 1, projection='polar')
for i, r_init in enumerate(r_init_ranges):
    color = colormap(normalize(r_init))
    print(color)
    params['r_init_multiplier'] = r_init
    N_max = 10
    if i == 0:
        N_max = 9.9
    c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = N_max, gamma=1)
    c.run_background_eq_of_motion()
    phi = c.sol['y'][3]
    theta = c.sol['y'][4]
    omega_p = c.get_omega_phi()
    point = 0
    for i in range(len(omega_p)):
        if omega_p[i] > 0.7:
            point = i
            break

    p = params['p']
    req = np.power(p*params['alpha']**2 /(6*params['m']**2*(params['V0'] - params['alpha']*theta)), 1/(2+p))
    print(point, len(omega_p))
    ax.plot(theta/30, phi, color=color)
    ax.plot(theta[point]/30, phi[point], 'go', color=color, linewidth=5, markersize=14)
ax.plot(np.linspace(0, 2*pi, 100), np.ones(100)*params['r0'], label=r'$r_0$')
#plt.polar(np.linspace(0, 2*pi, 100), req*np.ones(100), label=r'$r_{eq}$', c='k', linestyle='dashed', linewidth=3)
ax.plot(theta/30, req, label=r'$r_{eq}$', color='black', linestyle='dashed', linewidth=3)

#trans, _ , _ = ax.get_xaxis_text1_transform(-10)
#ax.text(np.deg2rad(45), -0.0, r"$\theta$", transform=trans, 
#         rotation=45-90, ha="center", va="center")

ax.set_rgrids([0.0025,0.005], angle=300)

label_position=ax.get_rlabel_position()
ax.text(np.radians(label_position+45),ax.get_rmax()/2.,r'$r$',
        rotation=-10,ha='center',va='center')

trans, _ , _ = ax.get_xaxis_text1_transform(0)
print(trans)
ax.text(np.deg2rad(90), 0.05, r"$\theta$", transform=trans, 
         rotation=0, ha="center", va="center")

plt.legend(loc='upper right')
plt.grid(True)
plt.box(on=None)
plt.xticks([])
fig.set_tight_layout(True)
fig.set_size_inches(7, 4)
plt.show()