import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation, PillowWriter

# This code has been created by explain how damping vibration works with different natural frequency. To get further information I fervently recommend you take a look inside my github
# For the one hand, the spring equation when the damping is proporcional to velocity ma + gammav + Kx = 0, for the other hand, the pendulum equation is theta'' - gamma/(L^2 M)theta' - g/L sin(theta)
class Damping_vibration():
    
    g = 9.81

    def __init__(self, x0, y0, v0, m, gamma, t_max = 15, dt = 0.01, pendulum = False, anim = False):

        #Physics parameter
        self.x0 = x0 # Initial position
        self.y0 = y0 # Initial position
        self.v0 = v0 # Angular or linear velocity
        self.m = m
        self.gamma = gamma
        self.pendulum = pendulum
        self.anim = anim

        #Times
        self.t = 0 
        self.t_max = t_max
        self.dt = dt

        # Stored histories
        self.Ek, self.Ep, self.Em = [], [], []
        self.y_hist, self.t_hist = [], []

    def run(self, pendulum, gamma):
        
         if pendulum:
            self.L = input('Pendulum length: ')
            approx = input('')
            self.betta = gamma/(self.L**2 * self.m)
            self.omega0 = np.sqrt(self.g/self.L)
            self.x_hist = []
            self.x_hist, self.y_hist, self.t_hist = self.Pendulum(self.betta, self.omega0, approx = False) # 2-Dimensional pendulum
         else:
            self.K = input('Elastic constant: ')
            self.betta = gamma/(2 * self.m)
            self.omega0 = np.sqrt(self.K / self.m)
            self.y_hist, self.t_hist = self.Spring(self.betta, self.omega0) # 1-Dimensional spring

class Pendulum(Damping_vibration):

    def __init__(self, theta0, omega0, m, L, betta, approx, *kwargs):
        super().__init__(theta0, omega0, m, *kwargs)

        self.L = L
        self.approx = bool(approx)

        # To ensure that any variable has been defined

        if not hasattr(self, 't'): self.t = 0.0 

        if not hasattr(self, 'dt'): self.dt = 1e-2 

        if not hasattr(self, 't_max'): self.t_max = 15.0


        self.t_hist = [self.t]
        self.omega_hist = [self.omega]
        self.theta_hist = [self.theta]

    def step(self):

        self.omega, self.theta = [], []

        while self.t_max > self.dt:
            if self.approx:
                omega_new, theta_new = self.Euler(self.dt)
            else: 
                if getattr(self, 'betta', 0) > 30:
                    omega_new, theta_new = self.Crank_Nicolson(self.dt)
                else:
                    omega_new, theta_new = self.rk4(self.dt) 
            
            # Update current postion
            self.dt += self.dt
            self.t_hist.append(self.dt)
            self.omega_hist.append(omega_new)
            self.theta_hist.append(theta_new)

            # Asign current state
            self.theta = theta_new
            self.omega = omega_new

    def Euler(self, dt):

        if abs(self.theta) > np.deg2rad(15):
            raise ValueError("The angle cannot be exceed 15 degrees")
            
        alpha = - self.betta * self.omega -(self.g / self.L) * self.theta

        omega_new += alpha * self.dt
        theta_new += self.omega * self.dt

        return omega_new, omega_new
    
    def rk4(self, dt):

        def f_theta(theta, omega):
            return omega
        def f_omega(theta, omega):
            return - self.betta * omega - (self.g / self.L) * np.sin(theta)
            
        # k1
        k1_theta = f_theta(self.theta, self.omega)
        k1_omega = f_omega(self.theta, self.omega)

        # k2
        k2_theta = f_theta(self.theta + 0.5 * dt * k1_theta, self.omega + 0.5 * dt * k1_omega)
        k2_omega = f_omega(self.theta + 0.5 * dt * k1_theta, self.omega + 0.5 * dt * k1_omega)

        # k3
        k3_theta = f_theta(self.theta + 0.5 * dt * k2_theta, self.omega + 0.5 * dt * k2_omega)
        k3_omega = f_omega(self.theta + 0.5 * dt * k2_theta, self.omega + 0.5 * dt * k2_omega)

        # k4
        k4_theta = f_theta(self.theta + dt * k3_theta, self.omega + dt * k3_omega)
        k4_omega = f_omega(self.theta + dt * k3_theta, self.omega + dt * k3_omega)

        theta_new += (dt / 6) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta)
        omega_new += (dt / 6) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega)

        return theta_new, omega_new
    
    def Crank_Nicolson(self, dt):

        x_new = self.theta
        v_new = self.omega
        max_iter = 10
        tol = 1e-6

        for k in range(max_iter):

            F1 = x_new - self.theta - dt/2*(self.omega + v_new) 
            F2 = v_new - self.omega - dt/2*(-self.betta * (self.omega + v_new) - (self.g/self.L) * (np.sin(self.theta) + np.sin(x_new))) 

            # Jacobian
            dF1_dx = 1 
            dF1_dv = -dt/2 
            dF2_dx = -dt/2 * (-(self.g/self.L) * np.cos(x_new)) 
            dF2_dv = 1 - dt/2 * (-self.betta)
            
            # Newton-Raphson update
            det = dF1_dx * dF2_dv - dF1_dv * dF2_dx

            if abs(det) < 1e-16: 
                raise RuntimeError("Jacobian singular en Crank-Nicolson")
            
            dx = (dF2_dv * F1 - dF1_dv * F2) / det
            dv = (dF1_dx * F2 - dF2_dx * F1) / det
            
            x_new -= dx
            v_new -= dv
            
            if abs(dx) < tol and abs(dv) < tol:
                break


        return x_new, v_new        
    
class Spring(Damping_vibration):