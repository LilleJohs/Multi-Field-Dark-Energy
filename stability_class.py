import numpy as np
from numpy import sqrt, exp, cos, sin
import scipy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

class MultiFieldDarkEnergy():
    def __init__(self, N_min = 0, N_max = 2, params=None, gamma=1, potential='exp', metric='squared'):
        if params == None:
            self.params = {
                'V0': 2.137,
                'm': 30,
                'r0': 7*1e-4,
                'alpha': 2*1e-3,
                'a': 2,
                'b': 2,
                'x_p_init': 1e-6,
                'x_t_init': 1e-7,
                'y_1_init': 1e-6,
                'r_init_multiplier': 1,
                'p': 2,
                'f0': 1,
                'beta': 1
            }
        else:
            self.params = params

        self.N_min = N_min
        self.N_max = N_max
        self.potential = potential
        self.metric = metric

        self.gamma = gamma 

        if self.potential == 'spinning' or self.potential == 'exp_spinning':
            if metric == 'exp':
                phi = self.params['r_init']
            else:
                phi = self.params['r0'] * self.params['r_init_multiplier']
            theta = 0
            x_p = self.params['x_p_init']
            y_1 = self.params['y_1_init']
            x_t = self.params['x_t_init']
            self.y0 = [x_p, x_t, y_1, phi, theta]
        else:
            x_p = self.params['x_p_init']
            y_1 = self.params['y_1_init']
            x_t = self.params['x_t_init']
            self.y0 = [x_p, x_t, y_1, 0.1, 0.01]

    def get_V_and_diff(self, phi, theta):
        if self.potential == 'axion':
            a = self.params['a']
            b = self.params['b']
            V = (1-cos(phi/a)) + (1-cos(theta/b))
            V_phi = sin(phi/a)/a
            V_theta = sin(theta/b)/b
        elif self.potential == 'spinning':
            V = self.params['V0'] - self.params['alpha']*theta + 0.5*self.params['m']**2 * (phi-self.params['r0'])**2
            V_phi = self.params['m']**2 * (phi-self.params['r0'])
            V_theta = - self.params['alpha']
        elif self.potential == 'exp':
            a = self.params['a']
            b = self.params['b']
            V = exp(a*phi + b*theta)
            V_phi = a*V
            V_theta = b*V
        elif self.potential == 'exp_lin':
            a = self.params['a']
            b = self.params['b']
            V = exp(a*phi) + exp(b*theta)
            V_phi = a*exp(a*phi)
            V_theta = b*exp(b*theta)
        elif self.potential == 'exp_spinning':
            V = self.params['V0']*exp(-self.params['alpha']*theta) + 0.5*self.params['m']**2 * (phi - self.params['r0'])**2 + self.params['cosmo_constant']
            V_phi = self.params['m']**2 * (phi - self.params['r0'])
            V_theta = - self.params['alpha']*self.params['V0']*exp(-self.params['alpha']*theta)
        return V, V_phi, V_theta

    def f_and_diff(self, phi):
        if self.metric == 'squared':
            return phi**2, 2*phi
        elif self.metric == 'cos':
            return cos(phi), -sin(phi)
        elif self.metric == 'flat':
            return 1, 0
        elif self.metric == 'r_p':
            return np.abs(phi)**self.params['p'],  self.params['p'] *np.abs(phi)**(self.params['p']-1)
        elif self.metric == 'exp':
            return self.params['f0']*np.exp(self.params['beta']*phi), self.params['f0']*self.params['beta']*np.exp(self.params['beta']*phi)

    def get_field_ricci(self, phi):
        if self.metric == 'squared':
            return 0
        elif self.metric == 'exp':
            return - self.params['beta']**2 / 2
        elif self.metric == 'r_p':
            return - self.params['p']*(self.params['p'] - 2)/ (2 * phi**2)
        else:
            return ''

    def background(self, N, y):
        x_p, x_t, y_1, phi, theta = y

        f, f_p = self.f_and_diff(phi)
        gamma = self.gamma

        if self.potential == 'exp':
            k3 = self.params['b']/sqrt(f)
            k2 = self.params['a']
        else:
            V, V_phi, V_theta = self.get_V_and_diff(phi, theta)
            k3 = V_theta/(sqrt(f)*V)
            k2 = V_phi/V
        k1 = f_p/f

        prime_x_p = 3*x_p*(x_p**2 + x_t**2 - 1) + sqrt(3/2)*(k1 * x_t**2 - k2 * y_1**2) + 3/2*gamma*x_p*(1 - x_p**2 - x_t**2 - y_1**2)
        prime_x_t = 3*x_t*(x_p**2 + x_t**2 - 1 ) - sqrt(3/2)*(k1 * x_p * x_t + k3 * y_1**2) + 3/2*gamma*x_t*(1 - x_p**2-x_t**2 - y_1**2)
        prime_y_1 = sqrt(3/2) * y_1 * (k2*x_p + k3*x_t) + 3/2*gamma*y_1*(1 - x_p**2 - x_t**2 - y_1**2) + 3*y_1*(x_p**2 + x_t**2)
        prime_phi = sqrt(6)*x_p
        prime_theta = sqrt(6)*x_t/sqrt(f)

        dydt = [prime_x_p, prime_x_t, prime_y_1, prime_phi, prime_theta]

        return dydt

    def run_background_eq_of_motion(self):
        sol = solve_ivp(fun = self.background, t_span = (self.N_min, self.N_max), y0 = self.y0, rtol=1.e-11)
        self.sol = sol

    def get_de_sitter_bound(self):
        phi = self.sol['y'][3]
        theta = self.sol['y'][4]
        V, V_phi, V_theta = self.get_V_and_diff(phi, theta)
        f, _ = self.f_and_diff(phi)

        return sqrt(V_phi**2 + V_theta**2 / f)/V

    def get_eq_of_state(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        y_1 = self.sol['y'][2]
        phi_squared = x_p**2 + x_t**2
        return (phi_squared - y_1**2) / (phi_squared + y_1**2)

    def get_omega_phi(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        y_1 = self.sol['y'][2]
        return x_p**2 + x_t**2 + y_1**2

    def get_field_derivative(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        N = self.sol['t']
        field_derivative = sqrt(6) * sqrt(x_p**2 + x_t**2)
        return field_derivative, np.trapz(field_derivative, N)
    
    def get_H(self):
        y_1 = self.sol['y'][2]
        phi = self.sol['y'][3]
        theta = self.sol['y'][4]
        V, _, _ = self.get_V_and_diff(phi, theta)
        return sqrt(V)/(sqrt(3)*y_1)
    
    def get_M_eff_squared(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        y_1 = self.sol['y'][2]
        phi = self.sol['y'][3]
        theta = self.sol['y'][4]

        f, _ = self.f_and_diff(phi)
        V, _, _ = self.get_V_and_diff(phi, theta)
        H = sqrt(V)/(sqrt(3)*y_1)
        dot_theta = sqrt(6)*H*x_t/sqrt(f)
        dot_phi = sqrt(6)*H*x_p
        dot_phi_squared = dot_phi**2 + f * dot_theta**2

        Vnn, omega = self.get_turning_rate()
        ricci = self.get_field_ricci(phi)
        M_eff_squared = (Vnn - omega**2 + ricci * dot_phi_squared/2)
        cs_s = 1/(1+4*omega**2/M_eff_squared)
        return M_eff_squared, cs_s

    def get_turning_rate(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        y_1 = self.sol['y'][2]
        r = self.sol['y'][3]
        theta = self.sol['y'][4]
        N = self.sol['t']

        f, f_r = self.f_and_diff(r)
        V, V_r, V_theta = self.get_V_and_diff(r, theta)
        H = sqrt(V)/(sqrt(3)*y_1)
        dot_theta = sqrt(6)*H*x_t/sqrt(f)
        dot_r = sqrt(6)*H*x_p

        Omega_squared = 1/f*(dot_r * V_theta - f * dot_theta * V_r)**2 / (dot_r**2 + f * dot_theta**2)**2

        phi_dot_norm = sqrt(dot_r**2 + f * dot_theta**2)
        if phi_dot_norm[0] == 0: phi_dot_norm[0]=1
        # Tangent vector upper indices
        tau_up_a = np.array([dot_r, dot_theta]) / phi_dot_norm
        # Normal vector upper indices
        
        N_up_a = np.array([sqrt(f) * tau_up_a[1], - 1/sqrt(f) * tau_up_a[0]])

        Vnn = N_up_a[0]**2 * self.params['m']**2 + N_up_a[1]**2 * self.params['alpha']**2*self.params['V0']*exp(-self.params['alpha']*theta) - N_up_a[1]**2 * (-f_r/2) * V_r - 2*N_up_a[0]*N_up_a[1] * (f_r/(2*f)) * V_theta
        
        return Vnn, sqrt(Omega_squared)


    def x_y_phase_plot(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        y_1 = self.sol['y'][2]
        phi = self.sol['y'][3]
        theta = self.sol['y'][4]

        f, f_r = self.f_and_diff(phi)
        V, V_phi, V_theta = self.get_V_and_diff(phi, theta)
        k3 = V_theta/(sqrt(f)*V)
        k2 = V_phi/V
        k1 = f_r/f

        plt.figure()
        plt.plot(x_p, y_1, label=r'$x_r$')
        plt.plot(x_t, y_1, label=r'$x_{\theta}$')
        plt.plot(-sqrt(1/6)*k3, y_1, label=r'$-k_3/\sqrt{6}$')
        plt.plot(sqrt(k2/k1)*y_1, y_1, label=r'$\sqrt{k_2/k_1}y_1$')
        #plt.plot(np.sqrt(x_p**2 + x_t**2), y_1, label=r'$\sqrt{x_{r}^2 + x_{\theta}^2}$')
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y_1$')
        plt.legend()
        plt.title('Phase Space x-y')
        #plt.savefig('img/x_y_x_1e-6_exp_spinning_r_2.pdf')
        plt.show()

    def plot_swampland_bound(self):
        x_p = self.sol['y'][0]
        x_t = self.sol['y'][1]
        y_1 = self.sol['y'][2]
        phi = self.sol['y'][3]
        theta = self.sol['y'][4]
        N = self.sol['t']

        f, _ = self.f_and_diff(phi)
        V, _, _ = self.get_V_and_diff(phi, theta)
        H = sqrt(V)/(sqrt(3)*y_1)
        dot_theta = sqrt(6)*H*x_t/sqrt(f)
        dot_phi = sqrt(6)*H*x_p
        dot_phi_squared = dot_phi**2 + f * dot_theta**2
        
        _, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 3))
        axes[0].plot(N, x_p, label=r'$\dot{\phi}/\sqrt{6}HM_p$')
        axes[0].plot(N, x_t, label=r'$\sqrt{f}\dot{\theta}/\sqrt{6}HM_p$')
        axes[0].plot(N, y_1, label=r'$\sqrt{V}/\sqrt{3}HM_p$')
        axes[0].legend()
        
        field_derivative, delta_phi = self.get_field_derivative()
        axes[1].plot(N, field_derivative, label=r'dPhi/dN')
        plt.title(delta_phi)
        axes[1].legend()

        N_now = 0
        for i in range(len(N)):
            if x_p[i]**2 + x_t[i]**2 + y_1[i]**2 >= 0.7:
                N_now = N[i]
                break

        _, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 3))
        axes[0].plot(N, x_p**2 + x_t**2, label=r'Kinetic')
        axes[0].plot(N, y_1**2, label=r'Potential')
        axes[0].plot(N, 1- x_p**2 - x_t**2 - y_1**2, label=r'$\Omega_M$')
        axes[0].legend()

        #axes[1].plot(N, self.get_eq_of_state(), label=r'$w$')
        axes[1].plot(N, self.get_de_sitter_bound(), label=r'De Sitter Bound')
        axes[1].plot(N, np.ones(len(N)), label=r'Lower Bound')
        #axes[1].plot(N, H, label=r'H')
        axes[1].set_yscale('log')
        axes[1].legend()

        _, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 3))
        if self.metric == 'r_p':
            r_eq = (self.params['p'] * self.params['alpha'] ** 2 * self.params['V0'] * exp(-self.params['alpha']*theta) / (6*self.params['m']**2) ) ** (1/(self.params['p']+2))
            axes[0].plot(N, r_eq, label=r'$r_{eq}$')
        #r_i = (self.params['p']*self.params['alpha']**2*self.params['V0']*exp(-self.params['alpha']*theta)/(6*self.params['m']**2))**(1/(self.params['p']+2))
        axes[0].plot(N, phi, label=r'$r$')
        
        axes[0].plot(N, self.params['r0']*np.ones(len(N)),label=r'$r_0$')
        #axes[0].plot(N, r_i, label=r'$r_i$')
        #axes[0].plot(N, dot_phi, label=r'$\dot{\phi}$')
        #axes[0].set_yscale('log')
        axes[0].legend()

        axes[1].plot(N, theta, label=r'$\theta$')
        axes[1].plot(N, dot_theta, label=r'$\dot{\theta}$')
        axes[1].set_yscale('log')
        #fig.tight_layout()
        axes[1].legend()

        plt.figure()
        w = self.get_eq_of_state()
        #weq = -1 + self.params['m'] * self.params['alpha']/(sqrt(3*self.params['V0']))*exp(self.params['alpha']*theta/2)
        if self.metric == 'r_p':
            weq = -1 + 2/(1 + self.params['p']/2 + self.params['p'] * self.params['V0'] * exp(-self.params['alpha'] * theta) / (self.params['m']**2 * r_eq**2))
            plt.plot(N, weq, label='Analytic')
        plt.plot(N, w, label='Equation of State')
        plt.vlines(N_now, min(w), max(w), 'r')
        
        plt.legend()

        plt.figure()
        Vnn, omega = self.get_turning_rate()
        ricci = self.get_field_ricci(phi)
        M_eff_squared = (Vnn - omega**2 + ricci * dot_phi_squared/2)
        plt.plot(N, Vnn, label='V_nn')
        plt.plot(N, -omega**2, label='-Omega^2')
        plt.plot(N, ricci * dot_phi_squared/2, label=r'$\mathcal{R}\dot{\phi}^2/2$')
        plt.plot(N, M_eff_squared, label='M^2_eff/a^2')
        plt.vlines(N_now, -1000, 1500, 'r')
        plt.legend()

        plt.figure()
        plt.plot(N, 1/(1+4*omega**2/M_eff_squared), label='c_s^2')
        plt.legend()

        f, f_p = self.f_and_diff(phi)
        V, V_phi, V_theta = self.get_V_and_diff(phi, theta)
        k3 = V_theta/(sqrt(f)*V)
        k2 = V_phi/V
        k1 = f_p/f
        plt.figure()
        plt.plot(N, k1, label='k1')
        plt.plot(N, k2, label='k2')
        plt.plot(N, k3, label='k3')
        plt.legend()
        #plt.yscale('log')
        plt.show()

       