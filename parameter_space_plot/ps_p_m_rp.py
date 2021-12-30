import sys
sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np
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

white = '#ffffff'
base_colors = sns.cubehelix_palette(3, reverse = False).as_hex()
base_colors[2] = base_colors[1]
base_colors[1] = white

print(base_colors)
cmap = colors.ListedColormap(base_colors)

my_cmap = np.array(cmap(np.arange(cmap.N)))
shaded_cmap = np.copy(my_cmap)
shaded_cmap[:, 2] += 0.2
shaded_cmap[1, :] = np.array([0.8, 0.8, 1, 1])

combined = np.concatenate((my_cmap, shaded_cmap))
print(combined)
cmap = colors.ListedColormap(combined)

bounds = [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
norm = colors.BoundaryNorm(bounds, cmap.N)

params = {
    'V0': 2.186,
    'm': 50,
    'r0': 7*1e-4,
    'alpha': 2*1e-3,
    'x_p_init': 0.0,
    'x_t_init': 0.0,
    'y_1_init': 1e-5,
    'p': 2,
    'r_init_multiplier': 1,
    'potential': 'spinning',
    'metric': 'r_p'
}

length = 100
list_accepted = np.zeros((length, length))

p_range = np.linspace(1.6, 3.2, length)
m_range = np.linspace(0, 500, length)
print(p_range)
print(m_range)
meff_over_p = np.zeros(length)

for j, m in enumerate(m_range):
    print(j, 'Mass:', m)
    for i, p in enumerate(p_range):
        params['p'] = p
        params['m'] = m

        #list_accepted[i, j]= get_sol(params)

fig, axs = plt.subplots(2)
axs[0].set_xlabel(r'$m \; [H_0]$')
axs[0].set_ylabel(r'$p$')
axs[0].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(p_range), np.amax(p_range)], aspect='auto')
#axs[0].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(p_range), np.amax(p_range)], aspect='auto',cmap='RdGy')
#axs[0].fill_between(m_range, meff_over_p, y2 = p_range[length-1], color='blue', alpha=0.4)

#np.save('list_accepted_m_p_r_p.npy', list_accepted)

params['p'] = 2
alpha_range = np.linspace(0, 0.02, length)
m_range = np.linspace(0, 500, length)
meff_over_alpha = np.zeros(length)

for j, m in enumerate(m_range):
    print(j, 'Mass:', m)
    for i, alpha in enumerate(alpha_range):
        params['alpha'] = alpha
        params['m'] = m
        
        cur_param = get_sol(params)
        list_accepted[i, j] = cur_param
axs[1].set_xlabel(r'$m \; [H_0]$')
axs[1].set_ylabel(r'$\alpha \; [H_0^2]$')
axs[1].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto')
#axs[1].imshow(list_accepted, origin = 'lower', extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto',cmap='RdGy')
#axs[1].fill_between(m_range, meff_over_alpha, y2 = alpha_range[length-1], color='blue', alpha=0.4)

fig.set_size_inches(7, 7)
#plt.savefig('test_r_p_m_p_alpha.pdf', bbox_inches = 'tight')
np.save('list_accepted_m_alpha_r_p.npy', list_accepted)

plt.show()