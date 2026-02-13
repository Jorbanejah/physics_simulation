import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Spring():

    g = 9.81

    def __init__(self, k, m, A, delta, y0, x0, Anim = False, t_max = 15):
        if k <= 0:
            raise ValueError("k must be positive")
        if m <= 0:
            raise ValueError("m must be positive")

        self.x = x0
        self.y0 = y0
        self.k = k 
        self.m = m 
        self.A = A 
        self.delta = delta 
        self.omega = np.sqrt(k/m) 
        self.t_max = t_max 
        self.Anima = Anim
        self.t = np.linspace(0, t_max, 100)
        self.time = 0

        if self.Anima:

            self.Ek, self.Ep, self.Ep_s, self.Em = [], [], [], []
            self.y_hist = []
            self.t_hist = []
            self.Animation()

    def compute(self, t):

        y =  self.A * np.cos(self.omega * t + self.delta) + self.y0
        vy =  - self.A * self.omega * np.sin(self.omega * t + self.delta)

        Ek = 0.5 * self.m * vy**2
        Ep_s = 0.5 * self.k * (y - self.y0)**2
        Em = Ek + Ep_s

        return y, Ek, Em, Ep_s
    
    def update(self, frame):
        
        y, Ek, Em, Ep_s = self.compute(frame)

        xs, ys = self.spring_shape(y)

        self.Ek.append(Ek)
        self.Em.append(Em)
        self.Ep_s.append(Ep_s)

        self.y_hist.append(y)
        self.t_hist.append(frame)

        self.point.set_data([self.x], [y])
        self.spring_line.set_data(xs, ys)

        self.point1.set_data(self.t_hist, self.y_hist)
        self.line_Ek.set_data(self.t_hist, self.Ek)
        self.line_Em.set_data(self.t_hist, self.Em)
        self.line_Ep_s.set_data(self.t_hist, self.Ep_s)

        if frame >= self.t_max - self.dt:
            self.anim.event_source.stop()

        return self.point, self.spring_line, self.point1, self.line_Ek, self.line_Em, self.line_Ep_s,

        
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
    
    def Animation(self):

        self.fig, self.ax = plt.subplots(1, 3, figsize= (14, 6))

        self.ax[0].axvline(self.x, color='black', lw=0.5)
        self.spring_line, = self.ax[0].plot([], [], color='black', linewidth=2)
        self.point, = self.ax[0].plot([], [], marker='o', markersize=12, color='red')
        self.point1, = self.ax[1].plot([], [], linewidth=2, color='blue')
        self.line_Ek, = self.ax[2].plot([], [], linewidth=2, color='red', label='Ek')
        self.line_Em, = self.ax[2].plot([], [], linewidth = 2, color='green', label='Em')
        self.line_Ep_s, = self.ax[2].plot([], [], linewidth = 2, color = 'blue', label ='Ep_s')
        

        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].set_xlim(- np.abs(self.x - 10), np.abs(self.x + 10))
        self.ax[0].set_ylim( -np.abs(self.y0*2), np.abs( self.y0 * 0.5))
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")
        self.ax[0].set_title("Spring")

        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim( - np.abs(self.y0 + self.y0 * 1.2),  np.abs(self.y0 - self.y0 *1.2))
        self.ax[1].axhline(self.y0, color='black', lw=0.5)
        self.ax[1].set_ylabel("y (m)")
        self.ax[1].set_xlabel("Time (s)")

        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].set_ylim(-np.abs(self.m * self.k * self.y0), np.abs(self.m * self.k * self.y0))
        self.ax[2].set_ylabel("E (J)")
        self.ax[2].set_xlabel("Time (s)")
        self.ax[2].legend()


        self.fig.tight_layout()

        self.dt = 0.1
        self.anim = FuncAnimation(self.fig, self.update, frames = np.arange(0, self.t_max, self.dt), interval = 50, blit = False)
        plt.show()

if __name__ == '__main__':
    animation = Spring(5, 1, 1,0, -2, 0, True)