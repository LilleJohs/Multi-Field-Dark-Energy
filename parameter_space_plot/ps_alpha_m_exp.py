import sys
sys.path.append("..")

from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from matplotlib import colors
import seaborn as sns
from get_sol_value import get_sol

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

base_colors = sns.cubehelix_palette(3, reverse = False).as_hex()
white = '#ffffff'
base_colors = sns.cubehelix_palette(3, reverse = False).as_hex()
base_colors[2] = base_colors[1]
base_colors[1] = white
cmap = colors.ListedColormap(base_colors)

my_cmap = np.array(cmap(np.arange(cmap.N)))
shaded_cmap = np.copy(my_cmap)
shaded_cmap[:, 2] += 0.2
shaded_cmap[1, :] = np.array([0.8, 0.8, 1, 1])

combined = np.concatenate((my_cmap, shaded_cmap))

cmap = colors.ListedColormap(combined)

bounds = [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
norm = colors.BoundaryNorm(bounds, cmap.N)

params = {
    'V0': 2.186,
    'm': 50,
    'r0': 0,
    'alpha': 1.5,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'beta': 800,
    'r_init': 0,
    'potential': 'spinning',
    'metric': 'exp'
}

length = 100
list_accepted = np.zeros((length, length))

beta_range = np.linspace(0, 3000, length)
m_range = np.linspace(0, 500, length)

meff_over_beta = np.zeros(length)
for j, m in enumerate(m_range):
    print(j, 'Mass:', m)
    for i, beta in enumerate(beta_range):
        params['beta'] = beta
        params['m'] = m
        
        list_accepted[i, j]= get_sol(params)

params['beta'] = 600
fig, axs = plt.subplots(2)
axs[0].set_xlabel(r'$m \; [H_0]$')
axs[0].set_ylabel(r'$\beta$')

axs[0].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(beta_range), np.amax(beta_range)], aspect='auto')
#axs[0].fill_between(m_range, meff_over_beta, y2 = beta_range[length-1], color='blue', alpha=0.4)

np.save('list_accepted_m_beta_exp.npy', list_accepted)
'''
alpha_range = np.linspace(0, 20, length)
m_range = np.linspace(0, 200, length)

for j, m in enumerate(m_range):
    print(j, 'Mass:', m)
    for i, alpha in enumerate(alpha_range):     
        params['alpha'] = alpha
        params['m'] = m
        list_accepted[i, j] = get_sol(params)

axs[1].set_xlabel(r'$m \; [H_0]$')
axs[1].set_ylabel(r'$\alpha \; [H_0^2]$')
axs[1].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto')
#axs[1].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto',cmap='RdGy')
#axs[1].fill_between(m_range, meff_over_alpha, y2 = alpha_range[length-1], color='blue', alpha=0.4)

fig.set_size_inches(7, 7)
#plt.savefig('test_alpha_m_exp.pdf', bbox_inches = 'tight')
np.save('list_accepted_m_alpha_exp.npy', list_accepted)
'''
plt.show()