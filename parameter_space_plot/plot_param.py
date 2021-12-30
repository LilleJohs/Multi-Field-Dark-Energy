import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import seaborn as sns

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
cmap = colors.ListedColormap(base_colors)

my_cmap = np.array(cmap(np.arange(cmap.N)))
shaded_cmap = np.copy(my_cmap)
shaded_cmap[:, 2] += 0.2
shaded_cmap[1, :] = np.array([0.8, 0.8, 1, 1])

combined = np.concatenate((my_cmap, shaded_cmap))

cmap = colors.ListedColormap(combined)

bounds = [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
norm = colors.BoundaryNorm(bounds, cmap.N)

length = 100
list_accepted = np.load('list_accepted_m_beta_exp.npy')
beta_range = np.linspace(0, 3000, length)
m_range = np.linspace(0, 500, length)

fig, axs = plt.subplots(2)
axs[0].set_title(r'$\alpha = 1.5 H_0^2$')
axs[0].set_xlabel(r'$m \; [H_0]$')
axs[0].set_ylabel(r'$\beta$')
axs[0].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(beta_range), np.amax(beta_range)], aspect='auto')

beta_scatter = [600, 1500]
m_scatter = [50, 400]
axs[0].scatter(m_scatter, beta_scatter, linewidths = 2.5, color='black')
axs[0].text(m_scatter[0] + 10,beta_scatter[0] + 50,'C')
axs[0].text(m_scatter[1] + 10,beta_scatter[1] + 50,'D')

alpha_range = np.linspace(0, 20, length)
m_range = np.linspace(0, 200, length)
list_accepted = np.load('list_accepted_m_alpha_exp.npy')

axs[1].set_title(r'$\beta = 600$')
axs[1].set_xlabel(r'$m \; [H_0]$')
axs[1].set_ylabel(r'$\alpha \; [H_0^2]$')
axs[1].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto')


fig.set_size_inches(7, 7)


p_scatter = [1.8, 2.3]
m_scatter = [400, 50]

fig, axs = plt.subplots(2)
p_range = np.linspace(1.6, 3.2, length)
m_range = np.linspace(0, 500, length)
list_accepted = np.load('list_accepted_m_p_r_p.npy')
axs[0].set_title(r'$\alpha = 2\cdot 10^{-3} H_0^2$')
axs[0].set_xlabel(r'$m \; [H_0]$')
axs[0].set_ylabel(r'$p$')
axs[0].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(p_range), np.amax(p_range)], aspect='auto')
axs[0].scatter(m_scatter, p_scatter, linewidths = 2.5, color='black')
axs[0].text(m_scatter[0] + 10,p_scatter[0]-0.02,'B')
axs[0].text(m_scatter[1] + 10,p_scatter[1]-0.05,'A')

alpha_range = np.linspace(0, 0.02, length)
m_range = np.linspace(0, 500, length)
list_accepted = np.load('list_accepted_m_alpha_r_p.npy')
axs[1].set_title(r'$p = 2$')
axs[1].set_xlabel(r'$m \; [H_0]$')
axs[1].set_ylabel(r'$\alpha \; [H_0^2]$')
axs[1].imshow(list_accepted, origin = 'lower', interpolation='nearest', cmap=cmap, norm=norm, extent=[np.amin(m_range), np.amax(m_range), np.amin(alpha_range), np.amax(alpha_range)], aspect='auto')

fig.set_size_inches(7, 7)
plt.show()