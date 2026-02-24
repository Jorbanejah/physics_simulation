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
        self.natural_frequency = np.sqrt(self.g / L)

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
        #if abs(self.theta) > np.deg2rad(15):
        #    raise ValueError("The angle cannot be exceed 15 degrees")

        alpha = -(self.g / self.L) * self.theta
        self.omega += alpha * self.dt
        self.theta += self.omega * self.dt
        return self.theta, self.omega

    def step_rk4(self):
        """Runge-Kutta 4"""

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
        E_max = 0.5 * self.m * (self.L)**2
        self.ax[2].set_ylim(-E_max * 1.1, E_max * 1.1)
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

class Spring:

    def __init__(self, k, m, A, delta, y0, x0, Anim = False, t_max=15, dt=0.05):

        if k <= 0:
            raise ValueError("k must be positive")
        if m <= 0:
            raise ValueError("m must be positive")

        # PHYSICAL PARAMETERS
        self.A = A
        self.delta = delta
        self.k = k
        self.m = m
        self.y0 = y0
        self.x = x0

        # Time
        self.t = 0
        self.dt = dt
        self.t_max = t_max

        # Natural frequency
        self.omega = np.sqrt(k/m)
        self.v_max = self.A * self.omega

        # Histories
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

    def compute_motion(self, t):

        y =  self.A * np.cos(self.omega * t + self.delta) + self.y0
        vy =  - self.A * self.omega * np.sin(self.omega * t + self.delta)

        return y, vy

    # ---------------------------------------------------------
    # ------------------  ENERGY  ---------------------------
    # ---------------------------------------------------------

    def compute_energy(self, y, v):

        Ek = 0.5 * self.m * v**2
        Ep = 0.5 * self.k * (y - self.y0)**2
        Em = Ek + Ep

        return Ek, Ep, Em

    # ---------------------------------------------------------
    # ------------------  SIMULATION WITHOUT ANIMATION ------------
    # ---------------------------------------------------------

    def run(self):

        while self.t <= self.t_max:

            y, v = self.compute_motion(self.t)
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
        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].set_xlim(- np.abs(self.x - 10), np.abs(self.x + 10))
        self.ax[0].set_ylim( -np.abs(self.y0*2), np.abs(self.y0 * 0.5))
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")
    
        
        self.point1, = self.ax[1].plot([], [], 'b-')
        self.ax[1].axhline(self.y0, color='black', lw=0.5)
        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_ylim( self.y0 - self.A - 2, self.y0 + self.A + 2)
        self.ax[1].set_ylabel("y (m)")
        self.ax[1].set_xlabel("Time (s)")


        self.line_Ek, = self.ax[2].plot([], [], 'r-', label='Ek')
        self.line_Ep, = self.ax[2].plot([], [], 'b-', label='Ep')
        self.line_Em, = self.ax[2].plot([], [], 'g-', label='Em')
        self.ax[2].axhline(0, color='black', lw=0.5)
        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].set_ylim(-np.abs(0.5 * self.m * self.v_max**2) - 2, np.abs(0.5 * self.m * self.v_max**2) + 2)
        self.ax[2].set_ylabel("E (J)")
        self.ax[2].set_xlabel("Time (s)")
        self.ax[2].legend()

        # Titles
        self.ax[0].set_title("Spring")
        self.ax[1].set_title("y(t)")
        self.ax[2].set_title("Energy")  

        self.fig.tight_layout()

        # Create the animation
        self.anim = FuncAnimation(
            self.fig, self._update_animation,
            frames=np.arange(0, self.t_max, self.dt),
            interval=50, blit=False
        )
        plt.show()

    def spring_shape(self, y_mass, n_coils=12, amplitude=0.3):
   
        y_top = 0
        ys = np.linspace(y_top, y_mass, n_coils*2)
        xs = []

        for i in range(len(ys)):
            if i == 0 or i == len(ys)-1:
                xs.append(self.x)
            else:
                xs.append(amplitude * (-1)**i + self.x)

        return xs, ys
    
    def _update_animation(self, frame):

        y, v = self.compute_motion(frame)
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

        if frame >= self.t_max - self.dt:
            self.anim.event_source.stop()

        return self.point, self.spring_line

if __name__ == '__main__':
    animation = Pendulum(L = 1, theta0=np.deg2rad(10), omega0=0, m=1, delta=0, t_max=10, dt=0.1, approx=False, animate=True)
    #animation = Spring(k = 2, m = 1, A = 1, delta = np.pi/2, y0 = -2, x0 = 1, Anim = True, t_max=15, dt=0.05)
