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
        self.x0 = x0
        self.y0 = y0
        self.v0 = v0
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

    
        if anim:
            self.setup_animation()
        else:
            self.run(pendulum, gamma)

    def run(self, pendulum, gamma):
        
         if pendulum:
            self.L = input('Pendulum length: ')
            approx = input()
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
class Spring(Damping_vibration):




