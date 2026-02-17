import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Pendulum:

    g = 9.81

    def __init__(self, L, theta0, omega0=0, m=1, delta=0,
                 t_max=15, dt=0.05, approx=False, animate=False):

        # Physics parameters
        self.L = L
        self.m = m
        self.theta = theta0          # angular position
        self.omega = omega0          # angular velocity
        self.delta = delta
        self.approx = approx

        # Time
        self.t = 0
        self.dt = dt
        self.t_max = t_max

        # Natural frequency
        self.omega0 = np.sqrt(self.g / L)

        self.t_hist = []
        self.x_hist = []
        self.Ek_hist = []
        self.Ep_hist = []
        self.Em_hist = []

        # Animation
        self.animate_flag = animate
        if animate:
            self._setup_animation()
        else:
            self.run()

    # ---------------------------------------------------------
    # ------------------  Pendulum Physics  -------------------
    # ---------------------------------------------------------

    def step_approx(self):
        """Approximation."""

        if abs(self.theta) > np.deg2rad(15):
            raise ValueError("The angle cannot be exceed 15 degrees")

        theta = self.theta * np.cos(self.omega0 * self.t + self.delta)
        omega = -self.theta * self.omega0 * np.sin(self.omega0 * self.t + self.delta)

        return theta, omega

    def step_rk4(self):
        """Runge–Kutta 4"""

        def f_theta(theta, omega):
            return omega

        def f_omega(theta, omega):
            return - (self.g / self.L) * np.sin(theta)

        dt = self.dt

        # k1
        k1_theta = f_theta(self.theta, self.omega)
        k1_omega = f_omega(self.theta, self.omega)

        # k2
        k2_theta = f_theta(self.theta + 0.5 * dt * k1_theta,
                       self.omega + 0.5 * dt * k1_omega)
        k2_omega = f_omega(self.theta + 0.5 * dt * k1_theta,
                       self.omega + 0.5 * dt * k1_omega)

        # k3
        k3_theta = f_theta(self.theta + 0.5 * dt * k2_theta,
                       self.omega + 0.5 * dt * k2_omega)
        k3_omega = f_omega(self.theta + 0.5 * dt * k2_theta,
                       self.omega + 0.5 * dt * k2_omega)

        # k4
        k4_theta = f_theta(self.theta + dt * k3_theta,
                       self.omega + dt * k3_omega)
        k4_omega = f_omega(self.theta + dt * k3_theta,
                       self.omega + dt * k3_omega)

        self.theta += (dt / 6) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta)
        self.omega += (dt / 6) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega)

        return self.theta, self.omega

    # ---------------------------------------------------------
    # ------------------  Energies  ---------------------------
    # ---------------------------------------------------------

    def compute_energy(self, theta, omega):
        Ek = 0.5 * self.m * (self.L * omega)**2
        Ep = self.m * self.g * self.L * (1 - np.cos(theta))
        Em = Ek + Ep
        return Ek, Ep, Em

    def run(self): 
        """Run the simulation without animation."""
        while self.t <= self.t_max:
            if self.approx:
                theta, omega = self.step_approx()
            else:
                theta, omega = self.step_rk4()
            Ek, Ep, Em = self.compute_energy(theta, omega)
            # Save history
            self.t_hist.append(self.t)
            self.x_hist.append(self.L * np.sin(theta))
            self.Ek_hist.append(Ek)
            self.Ep_hist.append(Ep)
            self.Em_hist.append(Em)
            self.t += self.dt
        
    # ---------------------------------------------------------
    # ------------------  ANIMATION  --------------------------
    # ---------------------------------------------------------

    def _setup_animation(self):
        self.fig, self.ax = plt.subplots(1, 3, figsize=(14, 6))

        # Pendulum
        self.line, = self.ax[0].plot([], [], 'k-', lw=2)
        self.point, = self.ax[0].plot([], [], 'ro', markersize=10)
        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].set_xlim(- self.L, self.L)
        self.ax[0].set_ylim(-self.L -1, 1)
        self.ax[0].set_title("Pendulum")
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")
        self.ax[0].set_aspect("equal", "box")

        # Position
        self.theta_line, = self.ax[1].plot([], [], 'b-')
        self.ax[1].axhline(0, color='black', lw=0.5)
        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim(-self.L, self.L)
        self.ax[1].set_title("x(t)")
        self.ax[1].set_xlabel("Time (s)")
        self.ax[1].set_ylabel("x (m)")

        # Energías
        self.Ek_line, = self.ax[2].plot([], [], 'r-', label="Ek")
        self.Ep_line, = self.ax[2].plot([], [], 'b-', label="Ep")
        self.Em_line, = self.ax[2].plot([], [], 'g-', label="Em")
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].set_ylim(-self.g, self.g)
        self.ax[2].set_title("Energy")
        self.ax[2].set_xlabel("Time (s)")
        self.ax[2].set_ylabel("E (J)")
        self.ax[2].legend()
        
        self.fig.tight_layout()

        self.anim = FuncAnimation(
            self.fig, self._update_animation,
            frames=np.arange(0, self.t_max, self.dt),
            interval=30, blit=False
        )
        plt.show()

    def _update_animation(self, frame):

        if self.approx:
            theta, omega = self.step_approx()
        else:
            theta, omega = self.step_rk4()

        # Energy
        Ek, Ep, Em = self.compute_energy(theta, omega)

        self.t_hist.append(self.t)
        self.x_hist.append(self.L * np.sin(theta))
        self.Ek_hist.append(Ek)
        self.Ep_hist.append(Ep)
        self.Em_hist.append(Em)

        # Cartesian coordinates
        x = self.L * np.sin(theta)
        y = -self.L * np.cos(theta)

        # Update pendulum position
        self.line.set_data([0, x], [0, y])
        self.point.set_data([x], [y])

        # Update plots
        self.theta_line.set_data(self.t_hist, self.x_hist)
        self.Ek_line.set_data(self.t_hist, self.Ek_hist)
        self.Ep_line.set_data(self.t_hist, self.Ep_hist)
        self.Em_line.set_data(self.t_hist, self.Em_hist)


        self.t += self.dt

        return self.line, self.point

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