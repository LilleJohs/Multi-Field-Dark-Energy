import sys
sys.path.append("..")
from stability_class import MultiFieldDarkEnergy
from shapely.geometry import Point, Polygon
import numpy as np
import matplotlib.pyplot as plt

cur_time = Polygon([(-1.05, 0.65), (-0.95, 0.75), (-0.95, 0.65), (-1.05, 0.75)])

def get_sol(params):
    c = MultiFieldDarkEnergy(params=params, N_min = 0, N_max = 9, gamma=1)
    c.run_background_eq_of_motion()

    H = c.get_H()
    # Today is when H=H_0 (which is one in unites of H_0)
    today_i = np.absolute(H-1).argmin()
    size = today_i + 1
    w = c.get_eq_of_state()[:size]
    omega = c.get_omega_phi()[:size]

    cur_param = 0.5
    for k in range(size):
        point = Point(w[k], omega[k])
        if cur_time.contains(point):
            # This solution has once been in omega=0.7 w=-1
            cur_param = 0.5
            break
        elif k == size-1:
            # This solution has NEVER been in omega=0.7 w=-1
            cur_param = 1
    if min(c.get_de_sitter_bound()[:size]) < 0.5:
        cur_param = 0
    if c.get_omega_phi()[len(H)-1] < 0.9 and c.get_eq_of_state()[len(H)-1] < -0.8:
        print('We should not be here')
        quit()
        #cur_param=-1
    M_eff_squared, cs_s = c.get_M_eff_squared()
    M_eff_squared = np.array(M_eff_squared)[:size]
    cs_s = np.array(cs_s[:size])

    # Check if the pert equations are valid when cs_s < 0
    negative_cs_s_indices = np.argwhere(cs_s < -1e-8)[:, 0]
    N = c.sol['t']-c.sol['t'][size-1]
    a = np.exp(N)

    lower = (H[:size] * a[:size])**2
    upper = M_eff_squared * a[:size] ** 2 / cs_s
    whenever = np.array(np.where((upper/lower)[negative_cs_s_indices] > 50**2))
    '''print(whenever.size)
    plt.figure()
    plt.plot(N[:size], upper / lower, label='upper / lower')
    plt.plot(N[:size], cs_s)
    plt.title('p: {}, m: {}, blue shade {}'.format(params['p'], params['m'], negative_cs_s_indices.size > 0 and whenever.size > 0))
    plt.legend()
    plt.show()'''
    if negative_cs_s_indices.size > 0 and whenever.size > 0:
        cur_param += 1.5
    ##if len(negative_cs_s_indices) > 0:
    #    cur_param += 1.5
    
    return cur_param