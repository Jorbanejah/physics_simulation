import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Pendulum:
    g = 9.81

    def __init__(self, L, delta, A, m, t_max=15):
        
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


        self.fig, self.ax = plt.subplots(1, 3, figsize=(14, 6))

        self.line_pendulum, = self.ax[0].plot([], [], color='blue')   
        self.point, = self.ax[0].plot([], [], marker='o', markersize=12, color='red')
        self.point1, = self.ax[1].plot([], [], linewidth=2, color='blue')
        self.line_Ek, = self.ax[2].plot([], [], linewidth=2, color='red', label='Ek')
        self.line_Ep, = self.ax[2].plot([], [], linewidth=2, color='blue', label='Ep')
        self.line_Em, = self.ax[2].plot([], [], linewidth = 2, color='green', label='Em')

        self.ax[2].legend()

        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim(-L - 0.5, L + 0.5)
        self.ax[1].axhline(0, color='black', lw=0.5)
        self.ax[1].set_ylabel("x (m)")
        self.ax[1].set_xlabel("Time (s)")

        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].set_ylim(-self.m * self.g * self.L* 2, self.m * self.g * self.L / 2)
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

        self.v = np.sqrt(vx**2 + vy**2)

        #Energy
        Ek = 0.5 * self.m * self.v**2
        Ep = self.m * self.g * (y)
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
            self.anim.event_source.stop()

        return self.point, self.point1, self.line_Ek, self.line_Ep, self.line_Em, self.line_pendulum,

    def animation(self):
        self.anim = FuncAnimation(self.fig, self.update, frames=np.arange(0, self.t_max, self.dt), interval=50, blit=False)
        plt.show()

class Spring:

    g = 9.81

    def __init__(self, k, m, A, delta, x0, y0, t_max=15):
        
        self.y0 = y0
        self.x = x0
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
        self.y_hist = []

        self.fig, self.ax = plt.subplots(1, 3, figsize= (14, 6))

        self.ax[0].axvline(self.x, color='black', lw=0.5)
        self.spring_line, = self.ax[0].plot([], [], color='black', linewidth=2)
        self.point, = self.ax[0].plot([], [], marker='o', markersize=12, color='red')
        self.point1, = self.ax[1].plot([], [], linewidth=2, color='blue')
        self.line_Ek, = self.ax[2].plot([], [], linewidth=2, color='red', label='Ek')
        self.line_Ep, = self.ax[2].plot([], [], linewidth=2, color='blue', label='Ep')
        self.line_Em, = self.ax[2].plot([], [], linewidth = 2, color='green', label='Em')
        

        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim( - np.abs(self.y0 + self.y0 * 1.2),  np.abs(self.y0 - self.y0 *1.2))
        self.ax[1].axhline(self.y0, color='black', lw=0.5)
        self.ax[1].set_ylabel("y (m)")
        self.ax[1].set_xlabel("Time (s)")

        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].axhline(( - np.abs(self.y0 * self.m * self.g)), color = 'black', lw = 0.5)
        self.ax[2].set_ylim(-np.abs(self.y0 + 2 * self.m * self.g* self.y0), 4)
        self.ax[2].set_ylabel("E (J)")
        self.ax[2].set_xlabel("Time (s)")
        self.ax[2].legend()

        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].set_xlim(- np.abs(self.x - 10), np.abs(self.x + 10))
        self.ax[0].set_ylim( -np.abs(self.y0*2), np.abs( self.y0 * 0.5))
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")
        self.ax[0].set_title("Spring")

        self.fig.tight_layout()

    def spring_shape(self, y_mass, n_coils=12, amplitude=0.3):
   
        y_top = 0
        ys = np.linspace(y_top, y_mass, n_coils*2)
        xs = []

        for i in range(len(ys)):
            if i == 0 or i == len(ys)-1:
                xs.append(0)
            else:
                xs.append(amplitude * (-1)**i)

        return xs, ys

    def update(self, frame):

        self.y =  self.A * np.cos(self.omega * self.time + self.delta) + self.y0
        vy =  - self.A * self.omega * np.sin(self.omega * self.time + self.delta)


        self.Ek.append(0.5 * self.m * vy**2)
        self.Ep.append(self.m * self.g * self.y)
        self.Em.append(self.Ek[-1] + self.Ep[-1])

        self.y_hist.append(self.y)
        self.t_hist.append(self.time)

        self.point.set_data([self.x], [self.y])

        xs, ys = self.spring_shape(self.y)
        self.spring_line.set_data(xs, ys)

        self.point1.set_data(self.t_hist, self.y_hist)
        self.line_Ek.set_data(self.t_hist, self.Ek)
        self.line_Ep.set_data(self.t_hist, self.Ep)
        self.line_Em.set_data(self.t_hist, self.Em)

        self.time += self.dt
        if self.time >= self.t_max:
            self.anim.event_source.stop()

        return self.point, self.spring_line, self.point1, self.line_Ek, self.line_Ep, self.line_Em,

    def animation(self):

        self.anim = FuncAnimation(self.fig, self.update, frames=np.arange(0, self.t_max, self.dt), interval=50, blit=False)
        plt.show()


if __name__ == "__main__":
    k = 0.8
    delta = 0
    A = 1  # radianes
    m = 1
    x0 = 0
    spring = Spring(k, m, A, delta, x0, y0 = - 2)
    spring.animation()

    A = np.pi*5/60 #15
    L = 1
    pendulum = Pendulum(L, delta, A, m)
    pendulum.animation()