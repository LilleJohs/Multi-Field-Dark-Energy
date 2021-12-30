import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from matplotlib.patches import ConnectionPatch

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
    'r_init_multiplier': 1,
    'p': 2,
    'cosmo_constant': 0,
    'potential': 'spinning',
    'metric': 'r_p'
}

def make_arrows(w_list, omega_list, params):
    arrow_x = np.zeros((len(w_list), len(omega_list)))
    arrow_y = np.zeros((len(w_list), len(omega_list)))

    for i, w in enumerate(w_list):
        print(i)
        for j, omega in enumerate(omega_list):
            # y^2 = omega*(1-w)/2
            # x^2 = omega*(1+w)/2
            params['y_1_init'] = sqrt(omega*(1-w)/2)
            #params['x_t_init'] = sqrt(omega*(1+w)/2)#/sqrt(2)
            params['x_p_init'] = sqrt(omega*(1+w)/2)#/sqrt(2)
            c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 2, gamma=1)

            c.run_background_eq_of_motion()
            #c.x_y_phase_plot()
            w_cur = c.get_eq_of_state()
            omega_cur = c.get_omega_phi()
            len_N = len(c.sol['t'])-1
            arrow_length = 1#sqrt((w_cur[len_N] - w)**2 + (omega_cur[len_N] - omega)**2)        
            
            arrow_x[j, i] = (w_cur[len_N] - w)/arrow_length
            arrow_y[j, i] = (omega_cur[len_N] - omega)/arrow_length
    return arrow_x, arrow_y

num = 9

red_box_x = np.array([[-1.005, -0.925], [-1.01, -0.74]])
red_box_y = np.array([[0.925, 1.005] , [0.74, 1.01]])

w_low = [-0.93, -0.75]

m_range = [50, 30]
p_range = [2, 3]
fig, axs = plt.subplots(2, 2)
for i, p in enumerate(p_range):
    params['p'] = p
    params['m'] = m_range[i]
    r_eq = (p*params['alpha']**2 /(6*params['m']**2 * params['V0']))**(1/(p+2))
    w_eq = -1 + 2/(1 + p/2 + p*params['V0']/ (r_eq**2* params['m']**2))

    w_list = np.linspace(start=-0.999, stop=w_low[i], num=num)
    omega_list = np.linspace(start=np.abs(w_low[i]), stop=0.9999, num=num)
    arrow_x, arrow_y = make_arrows(w_list, omega_list, params)
    if i == 0:
        axs[0, i].set_ylabel(r'$\Omega_{\phi}$')
    axs[0, i].quiver(w_list, omega_list, arrow_x, arrow_y)
    axs[0, i].plot(w_eq, 1, 'go', color='blue')
    axs[0, i].set_title(r'$p={{{}}}$'.format(p))
    axs[0, i].set_xlim(red_box_x[i, :])
    axs[0, i].set_ylim(red_box_y[i, :])
    w_list = np.linspace(start=-0.99, stop=0.95, num=num)
    omega_list = np.linspace(start=0.01, stop=0.99, num=num)
    arrow_x, arrow_y = make_arrows(w_list, omega_list, params)
    axs[1, i].set_xlabel(r'$w_{\phi}$')
    if i == 0:
        axs[1, i].set_ylabel(r'$\Omega_{\phi}$')
    axs[1, i].quiver(w_list, omega_list, arrow_x, arrow_y)
    axs[1, i].fill_between((red_box_x[i, 0], red_box_x[i, 1]), red_box_y[i, 0], red_box_y[i, 1], facecolor='red', alpha=0.2)
    
    con1 = ConnectionPatch(xyA=(red_box_x[i, 0], red_box_y[i, 0]), coordsA= axs[0, i].transData, 
                       xyB=(red_box_x[i, 0], red_box_y[i, 0]), coordsB= axs[1, i].transData, color = 'red')
    fig.add_artist(con1)

    con2 = ConnectionPatch(xyA=(red_box_x[i, 1], red_box_y[i, 0]), coordsA= axs[0, i].transData, 
                       xyB=(red_box_x[i, 1], red_box_y[i, 0]), coordsB= axs[1, i].transData, color = 'red')
    fig.add_artist(con2)
#axs[1].plot(w_eq, 1, 'go', color='blue')

#plt.savefig('p_init_arrow.pdf', bbox_inches = 'tight')
fig.set_size_inches(7, 7)

plt.show()