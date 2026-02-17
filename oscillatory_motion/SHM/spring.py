import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Spring():

    def __init__(self, k, m, A, delta, y0, x0, Anim=False, t_max=15, dt=0.01):

        if k <= 0:
            raise ValueError("k must be positive")
        if m <= 0:
            raise ValueError("m must be positive")

        # PHYSICAL PARAMETERS
        self.A = A
        self.delta = delta
        self.omega = np.sqrt(k/m)
        self.k = k
        self.m = m
        self.y0 = y0
        self.x = x0

        # initial conditions
        self.y = y0 + A * np.cos(delta)
        self.v = -A * np.sin(delta) * np.sqrt(k/m)

        # Time
        self.t = 0
        self.dt = dt
        self.t_max = t_max

        self.t_hist = []
        self.y_hist = []
        self.Ek_hist = []
        self.Ep_hist = []
        self.Em_hist = []

        # Animation
        self.Anim = Anim
        if Anim:
            self._setup_animation()
        else:
            self.run()  

    # ---------------------------------------------------------
    # ------------------ Spring physics  ------------------
    # ---------------------------------------------------------

    def step_euler(self):
        self.t += self.dt
        a = - (self.k / self.m) * (self.y - self.y0)
        self.v += a * self.dt
        self.y += self.v * self.dt
        return self.y, self.v

    # ---------------------------------------------------------
    # ------------------  ENERGÍAS  ---------------------------
    # ---------------------------------------------------------

    def compute_energy(self, y, v):
        Ek = 0.5 * self.m * v**2
        Ep = 0.5 * self.k * (y - self.y0)**2
        Em = Ek + Ep
        return Ek, Ep, Em

    # ---------------------------------------------------------
    # ------------------  SIMULACIÓN SIN ANIMACIÓN ------------
    # ---------------------------------------------------------

    def run(self):

        while self.t <= self.t_max:

            y, v = self.step_euler()
            Ek, Ep, Em = self.compute_energy(y, v)

            self.t_hist.append(self.t)
            self.y_hist.append(y)
            self.Ek_hist.append(Ek)
            self.Ep_hist.append(Ep)
            self.Em_hist.append(Em)

            self.t += self.dt

    # ---------------------------------------------------------
    # ------------------  ANIMATION  --------------------------
    # ---------------------------------------------------------

    def _setup_animation(self):

        self.fig, self.ax = plt.subplots(1, 3, figsize=(14, 6))

        self.spring_line, = self.ax[0].plot([], [], 'k-', lw=2)
        self.point, = self.ax[0].plot([], [], 'ro', markersize=10)
        self.ax[0].set_xlim(-0.5, 0.5)
        self.ax[0].set_ylim(-1.5, 0.5)
        
        self.point1, = self.ax[1].plot([], [], 'b-')
        self.ax[1].axhline(0, color='black', lw=0.5)
        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim(self.y0 - self.A*1.5, self.y0 + self.A*1.5)

        self.line_Ek, = self.ax[2].plot([], [], 'r-', label='Ek')
        self.line_Em, = self.ax[2].plot([], [], 'g-', label='Em')
        self.line_Ep, = self.ax[2].plot([], [], 'b-', label='Ep')
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].set_ylim(0, 1.5 * 0.5 * self.k * self.A**2)
        self.ax[2].legend()


        # Titles
        self.ax[0].set_title("Spring")
        self.ax[1].set_title("y(t)")
        self.ax[2].set_title("Energy")

        self.dt_anim = 0.05
        self.anim = FuncAnimation(
            self.fig, self._update_animation,
            frames=np.arange(0, self.t_max, self.dt_anim),
            interval=30, blit=False
        )
        plt.show()

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
    
    def _update_animation(self, frame):

        y, v = self.step_euler()
        Ek, Ep, Em = self.compute_energy(y, v)

        self.t_hist.append(self.t)
        self.y_hist.append(y)
        self.Ek_hist.append(Ek)
        self.Ep_hist.append(Ep)
        self.Em_hist.append(Em)

        xs, ys = self.spring_shape(y)

        self.point.set_data([self.x], [y])
        self.spring_line.set_data(xs, ys)

        self.point1.set_data(self.t_hist, self.y_hist)
        self.line_Ek.set_data(self.t_hist, self.Ek_hist)
        self.line_Ep.set_data(self.t_hist, self.Ep_hist)
        self.line_Em.set_data(self.t_hist, self.Em_hist)

        self.t += self.dt

        return self.point, self.spring_line

if __name__ == '__main__':
    animation = Spring(k = 10, m = 1, A = 0.5, delta = 0, y0 = -2, x0 = 0, Anim=True, t_max=15)