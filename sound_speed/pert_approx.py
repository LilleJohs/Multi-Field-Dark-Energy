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
    'V0': 2.186,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 2*1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'r_init_multiplier': 1, #r_i = x * r_0
    'p': 2,
    'potential': 'spinning',
    'metric': 'r_p'
}

fig, axs = plt.subplots(2)

h = 0.7
H0_in_mpc_inv = 1/4550

p_ranges = [1.8, 2.3]
m_ranges = [400, 50]
colormap = ListedColormap(sns.cubehelix_palette(10, reverse = False).as_hex())
normalize = mcolors.Normalize(vmin=np.min(p_ranges)*0.95, vmax=np.max(p_ranges)*1.05)
for i, p in enumerate(p_ranges):
    color = colormap(normalize(p_ranges[i]))
    params['p'] = p
    params['m'] = m_ranges[i]
    c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()
    N = c.sol['t']-8
    Vnn, omega = c.get_turning_rate()
    M_eff_s_over_a_s, cs_s = c.get_M_eff_squared()
    H = c.get_H()
    r = c.sol['y'][3]
    #cs_approx = (2-p+params['r0']/r * (p-1)) / (2+p-params['r0']/r * (p+1))
    a = np.exp(N)
    M_eff_s = a**2 * M_eff_s_over_a_s 
    #axs[0].plot(N, cs_s, label=r'$p = {{{}}}$'.format(p), color=color)
    #axs[0].plot(N, cs_approx, '--', label=r'Approx $p = {{{}}}$'.format(p), color=color)
    upper_bound = H0_in_mpc_inv * h * np.sqrt(M_eff_s/ cs_s)
    lower_bound = H0_in_mpc_inv * h * H * a
    axs[0].plot(N, upper_bound, c=color)
    axs[0].plot(N, lower_bound, '--', c=color)
    axs[0].fill_between(N, upper_bound / sqrt(50), lower_bound * sqrt(50), where=upper_bound / sqrt(50) > lower_bound * sqrt(50), alpha=0.4, color=color)
axs[0].set_yscale('log')
axs[0].set_xlim([-3, 2])
axs[0].set_ylabel(r'$k \; [h\mathrm{Mpc}^{-1}]$')
axs[0].legend([r'$|M_{\textrm{eff}} /c_s|$', r'$H a$'])
axs[0].set_title(r'$f(r)=r^p$')
axs[0].hlines(0, -8, 2, color='k', alpha=0.8)
axs[0].grid()

s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(p_ranges)

halfdist = (p_ranges[1] - p_ranges[0])/2.0
boundaries = np.linspace(p_ranges[0] - halfdist, p_ranges[-1] + halfdist, len(p_ranges) + 1)

#divider = make_axes_locat
cbar = fig.colorbar(s_map, ax=axs[0], ticks = p_ranges, boundaries = boundaries, spacing='proportional')
cbarlabel = r'$p$'
cbar.set_label(cbarlabel, fontsize=20)



params = {
    'V0': 2.186,
    'm': 50,
    'r0': 0,
    'alpha': 1.5,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'beta': 600,
    'r_init': 0,
    'potential': 'spinning',
    'metric': 'exp'
}

beta_ranges = [600, 1500]
m_ranges = [50, 400]
for i, beta in enumerate(beta_ranges):
    color = colormap(normalize(p_ranges[i]))
    params['beta'] = beta
    params['m'] = m_ranges[i]
    c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 10, gamma=1)
    c.run_background_eq_of_motion()

    H = c.get_H()
    # Today is when H=H_0 (which is one in unites of H_0)
    today_i = np.absolute(H-1).argmin()
    print('Time today:', c.sol['t'][today_i])
    N = c.sol['t']-c.sol['t'][today_i]

    M_eff_s_over_a_s, cs_s = c.get_M_eff_squared()
    H = c.get_H()
    a = np.exp(N)
    M_eff_s = M_eff_s_over_a_s * a**2
    upper_bound = H0_in_mpc_inv * h * np.sqrt(M_eff_s/ cs_s) 
    lower_bound = H0_in_mpc_inv * h* H * a 
    axs[1].plot(N, upper_bound, c=color)
    axs[1].plot(N, lower_bound, '--', c=color)
    axs[1].fill_between(N, upper_bound / 7, lower_bound * 7, where=upper_bound / 7 > lower_bound * 7, alpha=0.4, color=color)
axs[1].set_xlim([-3, 2])
axs[1].set_yscale('log')
axs[1].set_xlabel(r'$N$')
axs[1].set_ylabel(r'$k \; [h\mathrm{Mpc}^{-1}]$')
axs[1].legend([r'$|M_{\textrm{eff}} /c_s|$', r'$H a$'])
axs[1].set_title(r'$f(r)=e^{\beta r}$')
axs[1].hlines(0, -8, 2, color='k', alpha=0.8)
axs[1].grid()
fig.set_size_inches(7, 7)

colormap = ListedColormap(sns.cubehelix_palette(10, reverse = False).as_hex())
normalize = mcolors.Normalize(vmin=300, vmax=1800)
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(beta_ranges)

halfdist = (beta_ranges[1] - beta_ranges[0])/2.0
boundaries = np.linspace(beta_ranges[0] - halfdist, beta_ranges[-1] + halfdist, len(beta_ranges) + 1)

#divider = make_axes_locat

cbar = fig.colorbar(s_map, ax=axs[1], ticks = beta_ranges, boundaries = boundaries, spacing='proportional')
cbarlabel = r'$\beta$'
cbar.set_label(cbarlabel, fontsize=20)


plt.show()