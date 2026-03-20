import numpy as np

'''''
The following code describes the motion of two different systems: a nonlinear system and a linear system. Using, for instance, the Laplace transform, we can obtain the analytical solution of the linear system.
For the nonlinear case, we introduce a new symplectic method called the Verlet algorithm. In this version, the animation will only be generated for one system: the linear one, modeled as a Pohl pendulum.
The code has the following structure:
-    Integrators: RK4, Crank_Nicolson, and Verlet algorithms. Each of them is used for both systems.

-   Analytical solution: A simple function that returns the analytical solution of the linear system.

-   Main code: Implemented with three classes:
    -   DrivenOscillation: contains all parameters, flags, and general configuration.
    -   Linear: represents the linear system.
    -   Nonlinear: represents the nonlinear system

'''''
##
# ----------- Integrators --------------
##
def rk4(f, dt, q0, dq0, t, omega_0, omega, alpha, beta):
    """
    One RK4 step for a system:
        dq/dt  = v
        dv/dt  = f(q, v, t)

    Parameters
    ----------
    f   : function  -> acceleration function a = f(q, v)
    dt  : float     -> time step
    q0  : float     -> current position
    dq0 : float     -> current velocity

    Returns
    -------
    q_new, dq_new : floats
        Updated position and velocity
    """
    # k-values for velocity (dq/dt = v)
    k1_q = dq0
    k1_v = f(q0, dq0, t, omega_0, omega, alpha, beta)

    k2_q = dq0 + 0.5 * dt * k1_v
    k2_v = f(q0 + 0.5 * dt * k1_q, dq0 + 0.5 * dt * k1_v, t, omega_0, omega, alpha, beta)

    k3_q = dq0 + 0.5 * dt * k2_v
    k3_v = f(q0 + 0.5 * dt * k2_q, dq0 + 0.5 * dt * k2_v, t, omega_0, omega, alpha, beta)

    k4_q = dq0 + dt * k3_v
    k4_v = f(q0 + dt * k3_q, dq0 + dt * k3_v, t, omega_0, omega, alpha, beta)

    # RK4 update
    q_new  = q0  + (dt/6) * (k1_q + 2*k2_q + 2*k3_q + k4_q)
    dq_new = dq0 + (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)

    return q_new, dq_new

def crank_nicolson(f, dt, q0, dq0, t, omega_0, omega, alpha, beta):
    """
    One Crank_Nicolson step for the system:
        dx/dt = v
        dv/dt = f(x, v)

    Parameters
    ----------
    f   : function  -> acceleration function a = f(x, v)
    dt  : float     -> time step
    q0  : float     -> current position
    dq0  : float     -> current velocity

    Returns
    -------
    x_new, v_new : floats
        Updated position and velocity
    """

    #Explicit Euler
    x_new = q0 + dt * dq0
    v_new = dq0 + dt * f(q0, dq0, t, omega_0, omega, alpha, beta)

    max_iter = 10
    tol = 1e-6

    for _ in range(max_iter):

        # Crank–Nicolson equations
        F1 = x_new - q0 - 0.5 * dt * (dq0 + v_new)
        F2 = v_new - dq0 - 0.5 * dt * (f(q0, dq0) + f(x_new, v_new))

        # Jacobian entries
        dF1_dx = 1.0
        dF1_dv = -0.5 * dt

        # Partial derivatives of f(x, v)
        # Approximated numerically (finite differences)
        eps = 1e-6
        df_dx = (f(x_new + eps, v_new, t, omega_0, omega, alpha, beta) - f(x_new - eps, v_new, t, omega_0, omega, alpha, beta)) / (2 * eps)
        df_dv = (f(x_new, v_new + eps, t, omega_0, omega, alpha, beta) - f(x_new, v_new - eps, t, omega_0, omega, alpha, beta)) / (2 * eps)

        dF2_dx = -0.5 * dt * df_dx
        dF2_dv = 1.0 - 0.5 * dt * df_dv

        # Determinant
        det = dF1_dx * dF2_dv - dF1_dv * dF2_dx

        if abs(det) < 1e-14:
            raise RuntimeError("Jacobian is singular in Crank-Nicolson")

        # Newton correction
        dx = ( dF2_dv * F1 - dF1_dv * F2) / det
        dv = (-dF2_dx * F1 + dF1_dx * F2) / det

        # Update
        x_new -= dx
        v_new -= dv

        # Convergence check
        if abs(dx) < tol and abs(dv) < tol:
            break

    return x_new, v_new

def verlet(f, dt, q0, dq0, t, omega_0, omega, alpha, beta):

    '''
    One Verlet step for the system:
        dx/dt = v
        dv/dt = f(x, v)

    Parameters
    ----------
    f   : function  -> acceleration function a = f(x, v)
    dt  : float     -> time step
    q0  : float     -> current position
    dq0  : float     -> current velocity

    Returns
    -------
    x_new, v_new : floats
        Updated position and velocity
    '''
    x = np.array([q0])
    v = np.array([dq0])

    a = f(x, v, t, omega_0, omega, alpha, beta)
    x_new = x + v * dt + 0.5 * a * dt**2
    a_new = f(x, v, t, omega_0, omega, alpha, beta)
    v_new = v + dt/2 * (a_new + a)

    return x_new, v_new

##
#------------Functions and solution ----------
##


def linear_cos(x, v, t, omega_0, omega, alpha, beta):

    return - 2* beta * v - omega_0 * x + alpha * np.cos(omega * t)

def linear_sin(x, v, t, omega_0, omega, alpha, beta):

    return - 2* beta * v - omega_0 * x + alpha * np.sin(omega * t)

def nonlinear_cos(x, v, t, omega_0, omega, alpha, beta):

    return -2 * beta * v - omega_0 * np.sin(x) - alpha * np.cos(omega * t)

def nonlinear_sin(x, v, t, omega_0, omega, alpha, beta):

    return -2 * beta * v - omega_0 * np.sin(x) - alpha * np.sin(omega * t)

def analitics(A0, A, beta, t, omega, delta, chi):
    return A0 * np.exp(beta * t) * np.sin(omega * t + chi) + A * np.sin(omega* t - delta)


##
# ---------- Main code ---------------
##

class Driven_oscillation():

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
        self.st = dt

        #Set others parameters kwargs: k -- spring :: L -- pendulum
        self.system = system.lower()
        self.params = kwargs

        # Histories
        self.t_hist = []
        self.q_hist = []
        self.dq_hist = []

        self.Ek_hist = []
        self.Ep_hist = []
        self.Em_hist = []

    def run(self):
        
        self.F_external = self.params.get('F_external')
        
        if self.F_external is None:
            self.F_external = 'cos'

        elif self.F_external != 'cos' or self.F_external != 'sin':
            raise ValueError("The external force must be periodical: cos or sin")
            

        if self.system == 'linear':
            
            self.k = self.params.get("k")
            
            if self.k is None:
                raise ValueError("The system requiere a elastic constant k")
        
           
            self.model = Linear(y0 = self.q0, v0 = self.dq0, m = self.m, gamma = self.gamma, k = self.k, F0 = self.F0, dt = self.dt, t_max = self.t_max, omega = self.omega, F_external = self.F_external)

        elif self.system == 'nonlinear':
            
            self.L = self.params.get("L")

            if self.L is None:
                raise ValueError("The system requiere the pendulum length L")

            if -2 * np.pi > self.q0 > 2 * np.pi: #The angle must be in radians. Condition to change 

                self.q0 = np.deg2rad(self.q0)

                return self.q0
            
            self.model = Nonlinear(theta0 = self.q0, omega0 = self.dq0, m = self.m, gamma = self.gamma, L = self.L, F0 = self.F0, dt = self.dt, t_max = self.t_max, omega = self.omega, F_external = self.F_external)

class Linear():

    def __init__(self, y0, v0, m, gamma, k, F0, dt, t_max, omega, F_external):

        #Initial condition
        self.y0 = y0
        self.v0 = v0 
        
        #Parameters
        self.omega = omega # External force
        self.gamma = gamma
        self.k = k
        self.m = m

        self.beta = gamma/(2 * m) #Damping parameter
        self.omega2 = k/m #Natural frequency
        self.alpha = F0/m #Driven Force

        self.delta =np.arctan(2 * self.beta * self.omega / (self.omega2 - self.omega^2)) #gap between natural and external frequency

        self.dt = dt
        self.t_max = t_max

        self.F_external = F_external

    def run(self):
        numerical_methods = {
            "rk4": rk4,
            "CN": crank_nicolson,
            "Verlet": verlet
        }

        # Choose force function
        force_map = {
            "cos": linear_cos,
            "sin": linear_sin
        }

        force = force_map[self.F_external]


        # Prepare history dictionary
        self.history = {
            name: {"t": [], "q": [], "v": [], "Ek": [], "Ep": [], "Wp": []}
            for name in numerical_methods
        }

        # Loop over numerical methods
        for name, method in numerical_methods.items():

            # Reset initial conditions for each method
            q = self.y0
            dq = self.v0
            t = 0
            self.history[name]["Wp"].append(self.gamma * dq **2 * self.dt )
            while t < self.t_max:

                # Store current state

                # Energy
                self.history[name]["Ek"].append(0.5 * self.m * dq**2)
                self.history[name]["Ep"].append(0.5 * self.k * q**2)    
                wp = self.history[name]["Wp"[-1]] + self.gamma * dq**2 * self.dt
                self.history[name]["Wp"].append(wp)

                #Position
                self.history[name]["t"].append(t)
                self.history[name]["q"].append(q)
                self.history[name]["v"].append(dq)

                # Perform one integration step
                q, dq = method(force, q, dq, t, self.dt, t, self.omega2, self.omega, self.alpha, self.beta)

                t += self.dt

        self.analitical = {
            "x": [],
            "v": [],
            "t": []
        }
        t = 0
        y = self.y0
        v = self.v0

        while t < self.t_max:
            
            self.analitical["x"].append(y)
            self.analitical["y"].append(v)
            self.analitical["t"].append(t)

            y, v = analitics(y, v, self.beta, self.omega, self.delta, self.alpha, self.chi)

            t += self.dt

class Nonlinear():
    #At this point we have three kind of omega: omega0 --> initial radial velocity, omega02 ---> natural frequency, omega ---> external frequency
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
        self.omega02 = self.g/L
        self.alpha = F0/(m * L**2)

        #Time
        self.dt = dt
        self.t_max = t_max

        #External force
        self.omega = omega
        self.F_external = F_external
    
    def run(self):

        #Prepare motion history
        numerical_methods = {
            "Rk4": rk4,
            "Crank-nicolson": crank_nicolson,
            "Verlet": verlet
        }
        self.history = {
            name: {"t": [], "theta": [], "angular": [], "Ek": [], "Ep": [], "Wp": []}
            for name in numerical_methods
        }

        #External Force
        force_map = {
            "cos": nonlinear_cos, 
            "sin": nonlinear_sin
        }
        force = force_map[self.F_external]

        for name, method in numerical_methods.item():
            
            t = 0  
            theta = self.theta0
            angular = self.omega0
            self.history[name]["Wp"].append(self.gamma * angular**2 *self.dt)

            while t < self.t_max:

                # Energy
                self.history[name]["Ek"].append(0.5 * self.m * angular**2)
                self.history[name]["Ep"].append(self.g * self.L * theta**2)    
                wp = self.history[name]["Wp"[-1]] + self.gamma * angular**2 * self.dt
                self.history[name]["Wp"].append(wp)

                #Position
                self.history[name]["t"].append(t)
                self.history[name]["theta"].append(theta)
                self.history[name]["angular"].append(angular)

                theta, angular = method(force, self.dt, theta, angular, t, self.omega02, self.omega, self.alpha, self.beta) 
