import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.colors as mcolors

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

params = {
    'V0': 2.15,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1, #r_i = x * r_0
    'p': 2,
    'cosmo_constant': 0,
    'potential': 'spinning',
    'metric': 'r_p'
}

fig, axs = plt.subplots(2)

H0_in_mpc_inv = 1/4550

p_ranges = [0, 1, 2]
colormap = ListedColormap(sns.cubehelix_palette(10, reverse = False).as_hex())
normalize = mcolors.Normalize(vmin=np.min(p_ranges)-1, vmax=np.max(p_ranges)+1)
for i, p in enumerate(p_ranges):
    color = colormap(normalize(p_ranges[i]))
    params['p'] = p
    c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()
    N = c.sol['t']-8
    Vnn, omega = c.get_turning_rate()
    M_eff_s, cs_s = c.get_M_eff_squared()
    H = c.get_H()
    r = c.sol['y'][3]
    #cs_approx = (2-p+params['r0']/r * (p-1)) / (2+p-params['r0']/r * (p+1))
    a = np.exp(N)
    #axs[0].plot(N, cs_s, label=r'$p = {{{}}}$'.format(p), color=color)
    #axs[0].plot(N, cs_approx, '--', label=r'Approx $p = {{{}}}$'.format(p), color=color)
    upper_bound = H0_in_mpc_inv * np.sqrt(M_eff_s/ cs_s) * a
    lower_bound = H0_in_mpc_inv * H * a
    axs[0].plot(N, upper_bound, label=r'$|M_{eff} a /c_s| \;\;$' + r'$p={{{}}}$'.format(p))
    axs[0].plot(N, lower_bound, label=r'$H a \;\; p={{{}}}$'.format(p))
    axs[0].fill_between(N, upper_bound / 7, lower_bound * 7, where=upper_bound / 7 > lower_bound * 7, alpha=0.2)
axs[0].set_yscale('log')
axs[0].set_xlim([-3, 2])
axs[0].set_ylabel(r'$k \; [\mathrm{Mpc}^{-1}]$')
axs[0].legend()
axs[0].set_title(r'$f(r)=r^p$')
axs[0].hlines(0, -8, 2, color='k', alpha=0.8)
axs[0].grid()

params = {
    'V0': 2.15,
    'm': 50,
    'r0': 0,
    'alpha': np.sqrt(2),
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'beta': 100,
    'r_init': 0,
    'cosmo_constant': 0,
    'f0': 1,
    'potential': 'spinning',
    'metric': 'exp'
}

beta_ranges = [500, 1000, 1500]
colormap = ListedColormap(sns.cubehelix_palette(10, reverse = False).as_hex())
normalize = mcolors.Normalize(vmin=np.min(beta_ranges)-500, vmax=np.max(beta_ranges)+500)
for i, beta in enumerate(beta_ranges):
    color = colormap(normalize(beta_ranges[i]))
    params['beta'] = beta
    c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()
    N = c.sol['t']-8
    br = c.sol['y'][3] * beta
    M_eff_s, cs_s = c.get_M_eff_squared()
    H = c.get_H()
    a = np.exp(N)
    #axs[1].plot(N, cs_s, label=r'$\beta = {{{}}}$'.format(beta), color=color)
    #axs[1].plot(N, (1-br)/(1+br), '--', label=r'Approx $\beta = {{{}}}$'.format(beta), color=color)
    upper_bound = H0_in_mpc_inv * np.sqrt(M_eff_s/ cs_s) * a
    lower_bound = H0_in_mpc_inv * H * a
    axs[1].plot(N, upper_bound, label=r'$|M_{eff} a /c_s| \;\;$' + r'$\beta={{{}}}$'.format(beta))
    axs[1].plot(N, lower_bound, label=r'$H a \;\; \beta={{{}}}$'.format(beta))
    axs[1].fill_between(N, upper_bound / 7, lower_bound * 7, where=upper_bound / 7 > lower_bound * 7, alpha=0.2)
axs[1].set_xlim([-3, 2])
axs[1].set_yscale('log')
axs[1].set_xlabel(r'$N$')
axs[1].set_ylabel(r'$k \; [\mathrm{Mpc}^{-1}]$')
axs[1].legend()
axs[1].set_title(r'$f(r)=e^{\beta r}$')
axs[1].hlines(0, -8, 2, color='k', alpha=0.8)
axs[1].grid()
fig.set_size_inches(7, 7)

plt.show()