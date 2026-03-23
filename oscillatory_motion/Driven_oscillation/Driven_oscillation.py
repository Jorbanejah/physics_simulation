import numpy as np

"""
The following code describes the motion of two different systems: a nonlinear system and a linear system. Using, for instance, the Laplace transform, we can obtain the analytical solution of the linear system.
For the nonlinear case, we introduce a new symplectic method called the Verlet algorithm. In this version, the animation will only be generated for one system: the linear one, modeled as a Pohl pendulum.
The code has the following structure:
-    Integrators: RK4, Crank_Nicolson, and Verlet algorithms. Each of them is used for both systems.

-   Analytical solution: A simple function that returns the analytical solution of the linear system.

-   Main code: Implemented with three classes:
    -   DrivenOscillation: contains all parameters, flags, and general configuration.
    -   Linear: represents the linear system.
    -   Nonlinear: represents the nonlinear system

"""

##
# ----------- Integrators --------------
##
def rk4(f, t, y, dt, params=None):
    """
    Generic RK4 integrator

    Parameters
    ----------
    f : function
        f(t, y, params) -> dy/dt
    t : float
    y : np.ndarray
    dt : float
    params : optional

    Returns
    -------
    y_new : np.ndarray
    """

    k1 = f(t, y, params)
    k2 = f(t + dt/2, y + dt/2 * k1, params)
    k3 = f(t + dt/2, y + dt/2 * k2, params)
    k4 = f(t + dt, y + dt * k3, params)

    return y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

def crank_nicolson(f, t, y, dt, params=None, tol=1e-6, max_iter=10):

    """
    Generic implicit Crank-Nicolson
    """

    y_new = y + dt * f(t, y, params)  # Euler guess

    for _ in range(max_iter):

        F = y_new - y - 0.5 * dt * (
            f(t, y, params) + f(t + dt, y_new, params)
        )

        # Numerical Jacobian
        eps = 1e-6
        n = len(y)
        J = np.zeros((n, n))

        for i in range(n):
            dy = np.zeros(n)
            dy[i] = eps

            F_plus = (
                (y_new + dy) - y
                - 0.5 * dt * (f(t, y, params) + f(t + dt, y_new + dy, params))
            )

            F_minus = (
                (y_new - dy) - y
                - 0.5 * dt * (f(t, y, params) + f(t + dt, y_new - dy, params))
            )

            J[:, i] = (F_plus - F_minus) / (2 * eps)

        delta = np.linalg.solve(J, F)
        y_new -= delta

        if np.linalg.norm(delta) < tol:
            break

    return y_new

def velocity_verlet(accel, t, y, dt, params=None):
    """
    y = [q, v]
    accel(t, q, v, params)
    """

    q, v = y

    a = accel(t, q, v, params)

    q_new = q + v*dt + 0.5*a*dt**2
    a_new = accel(t + dt, q_new, v, params)

    v_new = v + 0.5*dt*(a + a_new)

    return np.array([q_new, v_new])

##
#------------Functions and solution ----------
##

def linear_system(t, y, p):
    """
    Linear system first order form

    y = [x, v]
    dy/dt = [v, a]
    """
    x, v = y

    # External force
    if p["type"] == "cos":
        force = np.cos(p["omega"] * t)
    else:
        force = np.sin(p["omega"] * t)

    dxdt = v
    dvdt = -2 * p["beta"] * v - p["omega0"] * x + p["alpha"] * force

    return np.array([dxdt, dvdt])

def nonlinear_system(t, y, p):
    """
    Non linear system first order form
    y = [x, v]
    dy/dt = [v, a]
    """
    x, v = y

    # External force
    if p["type"] == "cos":
        force = np.cos(p["omega"] * t)
    else:
        force = np.sin(p["omega"] * t)

    dxdt = v
    dvdt = -2 * p["beta"] * v - p["omega0"] * np.sin(x) + p["alpha"] * force

    return np.array([dxdt, dvdt])

def analytics(A0, A, beta, omega0, t, omega, delta, chi):

    omega_d = np.sqrt(omega0 - beta**2)
    
    return (A0 * np.exp(-beta * t) * np.sin(omega_d * t + chi) + A * np.sin(omega * t - delta))

##
# ---------- Main code ---------------
##

class DrivenOscillation():

    def __init__(self, q0, dq0, m, gamma, F0, omega, t = 15, dt = 0.01,  system = 'Linear', **kwargs):
        
        if m <= 0:
            raise ValueError('m must be positive')
        
        #Initial condition
        self.q0 =q0
        self.dq0 = dq0

        #Set parameters
        self.m = m
        self.gamma = gamma
        self.F0 = F0
        self.omega = omega #External frequency

        #Set times
        self.t_max = t
        self.dt = dt

        #Set others parameters kwargs: k -- spring :: L -- pendulum
        self.system = system.lower()
        self.params = kwargs

        if self.model not in ('linear', 'nonlinear'):
            raise ValueError("System must be linear or nonlinear")

    def run(self):
        
        #External force selection
        self.F_external = self.params.get('F_external')
        
        if self.F_external not in ('cos', 'sin'):
            raise ValueError("The external force must be periodical: cos or sin")
        
        #Model selection
        if self.system == 'linear':
            
            self.k = self.params.get("k")
            
            if self.k is None:
                raise ValueError("The system requiere a elastic constant k")
        
            self.model = Linear(y0 = self.q0, v0 = self.dq0, m = self.m, gamma = self.gamma, k = self.k, F0 = self.F0, dt = self.dt, t_max = self.t_max, omega = self.omega, F_external = self.F_external)

        else:
            self.L = self.params.get("L")

            if self.L is None:
                raise ValueError("The system requiere the pendulum length L")

            if  abs(self.q0) > 2*np.pi: #The angle must be in radians. Condition to change 

                raise ValueError("The angle theta must be in radians")
            
            self.model = Nonlinear(theta0 = self.q0, omega0 = self.dq0, m = self.m, gamma = self.gamma, L = self.L, F0 = self.F0, dt = self.dt, t_max = self.t_max, omega = self.omega, F_external = self.F_external)

        self.model.run()

        return self.model

class Linear:
    #At this point we have three kind of omega: omega_sq ---> natural frequency, omega ---> external frequency
    def __init__(self, y0, v0, m, gamma, k, F0, dt, t_max, omega, F_external):

        # Initial conditions
        self.y0 = y0
        self.v0 = v0 
        
        # Parameters
        self.m = m
        self.gamma = gamma
        self.k = k
        self.omega = omega  # external frequency

        # Derived parameters
        self.beta = gamma / (2 * m)
        self.omega_sq = k / m 
        self.alpha = F0 / m

        # Phase (steady-state)
        self.delta = np.arctan((2 * self.beta * self.omega) / (self.omega_sq - self.omega**2))
        self.chi = 0 

        # Time
        self.dt = dt
        self.t_max = t_max

        # External force type
        self.F_external = F_external

    def run(self):

        n_steps = int(np.ceil(self.t_max / self.dt))

        # Parameters for system
        params = {
            "beta": self.beta,
            "omega0": self.omega_sq,
            "omega": self.omega,
            "alpha": self.alpha,
            "type": self.F_external
        }

        # --- Integrators ---
        numerical_methods = {
            "rk4": rk4,
            "CN": crank_nicolson,
            "Verlet": velocity_verlet
        }

        # --- History ---
        self.history = {
            name: {"t": [], "q": [], "v": [], "Ek": [], "Ep": [], "Wp": []}
            for name in numerical_methods
        }

        # --- Acceleration adapter for Verlet ---
        def accel(t, q, v, p):
            return linear_system(t, np.array([q, v]), p)[1]

        # --- Loop over methods ---
        for name, method in numerical_methods.items():

            y = np.array([self.y0, self.v0])
            wp = 0

            for i in range(n_steps):
                t = i * self.dt
                q, dq = y

                # Energies
                Ek = 0.5 * self.m * dq**2
                Ep = 0.5 * self.k * q**2
                wp += self.gamma * dq**2 * self.dt

                # Store
                self.history[name]["t"].append(t)
                self.history[name]["q"].append(q)
                self.history[name]["v"].append(dq)
                self.history[name]["Ek"].append(Ek)
                self.history[name]["Ep"].append(Ep)
                self.history[name]["Wp"].append(wp)

                # Step
                if name == "Verlet":
                    y = method(accel, t, y, self.dt, params)
                else:
                    y = method(linear_system, t, y, self.dt, params)

        # --- Analytical solution ---
        self.analytical = {"x": [], "t": []}

        for i in range(n_steps):
            t = i * self.dt

            x = analytics(A0 = self.y0, A = self.alpha, beta = self.beta, omega0 = self.omega_sq, t=t, omega=self.omega, delta=self.delta, chi=self.chi)

            self.analytical["x"].append(x)
            self.analytical["t"].append(t)

class Nonlinear():
    #At this point we have three kind of omega: omega0 --> initial radial velocity, omega_sq ---> natural frequency, omega ---> external frequency
    g = 9.81

    def __init__(self, theta0, omega0, m, gamma, L, F0, dt, t_max, omega, F_external):
        #Initial conditions
        self.theta0 = theta0
        self.omega0 = omega0

        #Parameters
        self.m = m 
        self.L = L
        self.gamma = gamma
        self.F0 = F0 

        #Pendulum parameter
        self.beta = gamma/(m * L**2)
        self.omega_sq = self.g/L
        self.alpha = F0/(m * L**2)

        #Time
        self.dt = dt
        self.t_max = t_max

        #External force
        self.omega = omega
        self.F_external = F_external
    
    def run(self):
        
        n_steps = int(np.ceil(self.t_max/self.dt))

        params = {
            "beta": self.beta,
            "omega0": self.omega_sq,
            "omega": self.omega,
            "alpha": self.alpha,
            "type": self.F_external
        }

        numerical_methods = {
            "rk4": rk4,
            "Crank-nicolson": crank_nicolson,
            "Verlet": velocity_verlet
        }

        self.history = {
            name: {"t": [], "q": [], "v": [], "Ek": [], "Ep": [], "Wp": []}
            for name in numerical_methods
        }

        #Prepare motion history
        def accel(t, q, v, p):
            return nonlinear_system(t, np.array([q, v]), p)[1]

        # --- Loop over methods ---
        for name, method in numerical_methods.items():

            y = np.array([self.theta0, self.omega0])
            wp = 0

            for i in range(n_steps):
                t = i * self.dt
                q, dq = y

                # Energies
                Ek = 0.5 * self.m * dq**2
                Ep = self.m * self.g * self.L* (1 - np.cos(q))
                wp += self.gamma * dq**2 * self.dt

                # Store
                self.history[name]["t"].append(t)
                self.history[name]["q"].append(q)
                self.history[name]["v"].append(dq)
                self.history[name]["Ek"].append(Ek)
                self.history[name]["Ep"].append(Ep)
                self.history[name]["Wp"].append(wp)

                # Step
                if name == "Verlet":
                    y = method(accel, t, y, self.dt, params)
                else:
                    y = method(nonlinear_system, t, y, self.dt, params)