import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Pendulum:
    g = 9.81

    def __init__(self, L, delta, A, m, t_max=10):
        
        self.m = m
        self.L = L
        self.delta = delta
        self.A = A  # (rad)
        self.omega = np.sqrt(self.g / L)
        self.t_max = t_max
        self.time = 0
        self.dt = 0.05

        #Energy and position history
        self.Em, self.Ek, self.Ep = [], [], []
        self.t_hist = []
        self.x_hist = []


        self.fig, self.ax = plt.subplots(1, 3, figsize=(12, 4))

        self.line_pendulum, = self.ax[0].plot([], [], color='blue')   
        self.point, = self.ax[0].plot([], [], marker='o', markersize=12, color='red')
        self.point1, = self.ax[1].plot([], [], linewidth=2, color='blue')
        self.line_Ek, = self.ax[2].plot([], [], linewidth=2, color='red', label='Ek')
        self.line_Ep, = self.ax[2].plot([], [], linewidth=2, color='blue', label='Ep')
        self.line_Em, = self.ax[2].plot([], [], linewidth = 2, color='green', label='Em')

        self.ax[2].legend()

        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim(-L, L)
        self.ax[1].axhline(0, color='black', lw=0.5)
        self.ax[1].set_ylabel("x (m)")
        self.ax[1].set_xlabel("Time (s)")

        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].set_ylim(-self.m * self.g * self.L / 2, self.m * self.g * self.L / 2)
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].set_ylabel("E (J)")
        self.ax[2].set_xlabel("Time (s)")

        self.ax[0].set_xlim(-self.L - 1, self.L + 1)
        self.ax[0].set_ylim(-self.L - 1, 1)
        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].axvline(0, color='black', lw=0.5)
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")
        self.ax[0].set_title("Pendulum")

        # prevent the overlapping of subplots
        self.fig.tight_layout()

    def update(self, frame):

        #Position
        theta = self.A * np.cos(self.omega * self.time + self.delta)
        dtheta = -self.A * self.omega * np.sin(self.omega * self.time + self.delta)

        x = self.L * np.sin(theta)
        y = - self.L * np.cos(theta)

        vx = self.L * dtheta * np.cos(theta)
        vy = self.L * dtheta * np.sin(theta)

        v = np.sqrt(vx**2 + vy**2)

        #Energy
        Ek = 0.5 * self.m * v**2
        Ep = self.m * self.g * (y + self.L)
        Em = Ek + Ep

        self.Ek.append(Ek)
        self.Ep.append(Ep)
        self.Em.append(Em)
        self.t_hist.append(self.time)
        self.x_hist.append(x)

        self.point.set_data([x], [y])
        self.line_pendulum.set_data([0, x], [0, y])
        self.point1.set_data([self.t_hist], [self.x_hist])
        self.line_Ek.set_data([self.t_hist], [self.Ek])
        self.line_Ep.set_data([self.t_hist], [self.Ep])
        self.line_Em.set_data([self.t_hist], [self.Em])

        self.time += self.dt
        
        if self.time >= self.t_max:
            raise StopIteration("Maximum time reached")

        return self.point, self.point1, self.line_Ek, self.line_Ep, self.line_Em, self.line_pendulum,

    def animation(self):
        self.anim = FuncAnimation(self.fig, self.update, frames=np.arange(0, self.t_max, self.dt), interval=50, blit=False)
        plt.show()

class spring:
    def __init__(self, k, m, A, delta, t_max=10, y):
        self.y = y
        self.k = k 
        self.m = m 
        self.A = A 
        self.delta = delta 
        self.omega = np.sqrt(k/m) 
        self.t_max = t_max 
        self.time = 0 
        self.dt = 0.05

        self.Em, self.Ek, self.Ep = [], [], []
        self.t_hist = []
        self.x_hist = []

        self.fig, self.ax = plt.subplots(1, 3, figsize= (10, 6))
        self.point, = self.ax[0].plot([], [], marker='o', markersize=12, color='red')
        self.point1, = self.ax[1].plot([], [], linewidth=2, color='blue')
        self.line_Ek, = self.ax[2].plot([], [], linewidth=2, color='red', label='Ek')
        self.line_Ep, = self.ax[2].plot([], [], linewidth=2, color='blue', label='Ep')
        self.line_Em, = self.ax[2].plot([], [], linewidth = 2, color='green', label='Em')

        self.ax[2].legend()

        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].axhline(0, color='black', lw=0.5)
        self.ax[1].set_ylabel("x (m)")
        self.ax[1].set_xlabel("Time (s)")

        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].set_ylabel("E (J)")
        self.ax[2].set_xlabel("Time (s)")

        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].axvline(y, color='black', lw=0.5)
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")
        self.ax[0].set_title("Pendulum")

    def update(self, frame):
        pass
    def animation(self, frame):

        self.anim = FuncAnimation(self.fig, self.update, frames=np.arange(0, self.t_max, self.dt), interval=50, blit=False)
        plt.show()


if __name__ == "__main__":
    L = 10
    delta = 0
    A = np.pi/2  # radianes
    m = 1
    pendulum = Pendulum(L, delta, A, m)
    pendulum.animation()
