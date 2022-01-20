import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')
plt.rc('font', size=15)
plt.rc('figure', autolayout=True)
plt.rc('axes', titlesize=16, labelsize=17)
plt.rc('lines', linewidth=2, markersize=6)
plt.rc('legend', fontsize=15)
plt.rc('text', usetex=True)

def solve_for_req(params, theta):
    func = lambda req : req**3 + 2*params['V0']/params['m']**2 * np.exp(-params['alpha'] * theta) * req - (params['V0']*params['alpha']*np.exp(-params['alpha']*theta)/(params['m']**2))**2 * params['beta']/3 * np.exp(-params['beta']*req)
    sol = fsolve(func, 0.01)
    return sol

params = [{
    'V0': 2.186,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 2*1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1,
    'p': 2,
    'potential': 'spinning',
    'metric': 'r_p'
}, {
    'V0': 2.186,
    'm': 50,
    'r0': 0,
    'alpha': 3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'beta': 600,
    'r_init': 0,
    'potential': 'spinning',
    'metric': 'exp'
}]

r_lim = [-5, 1]
title_l = [r'$f(r)=r^p$', r'$f(r)=e^{\beta r}$']

r = np.logspace(r_lim[0], r_lim[1], num=100)
theta_times_alpha_l = np.array([[0, 1, 2], [0, 1, 2]])
print(theta_times_alpha_l.shape)
fig, axs = plt.subplots(2)

for j in range(2):
    colormap = ListedColormap(sns.cubehelix_palette(len(theta_times_alpha_l[j, :]), reverse = False).as_hex())
    normalize = mcolors.Normalize(vmin=np.min(theta_times_alpha_l[j, :]), vmax=np.max(theta_times_alpha_l[j, :]))

    c = MultiFieldDarkEnergy(params=params[j], N_min = 0, N_max = 10, gamma=1)
    for i, theta_times_alpha in enumerate(theta_times_alpha_l[j, :]):
        color = colormap(normalize(theta_times_alpha_l[j, i]))
        theta = theta_times_alpha / params[j]['alpha']
        V, V_r, V_theta = c.get_V_and_diff(r, theta)
        f, _ = c.f_and_diff(r)
        de_sitter_c = sqrt(V_r**2 + f**(-1)*V_theta**2)/V
        if j==0:
            param = params[0]
            req = np.power(param['p']*param['alpha']**2 * (param['V0'] -param['alpha']*theta)/(6*param['m']**2), 1/(2+param['p']))
        else:
            req = solve_for_req(params[j], theta)
        axs[j].axvline(x=req, color=color, linestyle='--', linewidth=2)

        axs[j].plot(r, de_sitter_c, color=color)
        
    axs[j].axvline(x=params[j]['r0'], color='grey', linestyle='-.')
    axs[j].hlines(1, xmin=10**(r_lim[0]), xmax=10**(r_lim[1]), color='grey')

    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    s_map.set_array(theta_times_alpha_l[j,:])

    halfdist = (theta_times_alpha_l[j, 1] - theta_times_alpha_l[j, 0])/2.0
    boundaries = np.linspace(theta_times_alpha_l[j, 0] - halfdist, theta_times_alpha_l[j, -1] + halfdist, len(theta_times_alpha_l[j, :]) + 1)

    #divider = make_axes_locatable(axs[0])
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    #cbar=fig.colorbar(im1, cax=cax, orientation='vertical')

    cbar = fig.colorbar(s_map, ax=axs[j], ticks = theta_times_alpha_l[j, :], boundaries = boundaries, spacing='proportional')
    cbarlabel = r'$\alpha \theta$'
    cbar.set_label(cbarlabel, fontsize=20)

    axs[j].set_title(title_l[j])
    axs[j].set_xscale('log')
    axs[j].set_yscale('log')
    axs[j].set_xlabel(r'$r$')
    axs[j].set_ylabel(r'$|\nabla V| / V$')
    axs[j].set_xlim([10**(r_lim[0]), 10**(r_lim[1])])
fig.set_size_inches(7, 7)
plt.show()