'''''
DampedVibration
│
├── Pendulum
│     ├ solve()
│     ├ rk4()
│     └ crank_nicolson()
│
└── Spring
      ├ solve()
      └ rk4()       
'''''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# This code has been created by explain how damping vibration works with different natural frequency. To get further information I fervently recommend you take a look inside my github
# For the one hand, the spring equation when the damping is proporcional to velocity ma + gammav + Kx = 0, for the other hand, the pendulum equation is theta'' - gamma/(L^2 M)theta' - g/L sin(theta)

class DampedVibration:
    """
    General wrapper for damped mechanical systems.
    Supported systems:
        - Pendulum
        - Spring
    """
    def __init__(self, q0, dq0, m, gamma, t_max=15, dt=0.01, system="pendulum", animate=False, **kwargs):

        # Initial conditions: generalize coordenates
        self.q0 = q0
        self.dq0 = dq0

        # Physical parameters
        self.m = m
        self.gamma = gamma

        # Time parameters
        self.t_max = t_max
        self.dt = dt

        # Configuration
        self.system = system.lower()
        self.animate = animate

        # Extra parameters (L or k)
        self.params = kwargs

        # Histories
        self.t_hist = []
        self.q_hist = []
        self.dq_hist = []

        self.Ek_hist = []
        self.Ep_hist = []
        self.Em_hist = []

        self.model = None

    def run(self):

        if self.system == "pendulum":
            
            self.L = self.params.get("L")
            approx = self.params.get("approx", False)

            if self.L is None:
                raise ValueError("Pendulum requires length L")

            self.model = Pendulum(theta0=self.q0, omega0= self.dq0, m=self.m, gamma=self.gamma, L=self.L, dt=self.dt, t_max=self.t_max, approx=approx)

        elif self.system == "spring":

            self.k = self.params.get("k")

            if self.k is None:
                raise ValueError("Spring requires elastic constant k")

            self.model = Spring(y0=self.q0, v0= self.dq0, m=self.m, gamma=self.gamma, k=self.k, dt=self.dt, t_max=self.t_max)

        else:

            raise ValueError('The class only works with pendulum or spring system')
        
        # Run simulation
        self.model.solve()

        # unified interface
        self.t_hist = self.model.t_hist
        self.q_hist = self.model.q_hist
        self.dq_hist = self.model.dq_hist

        self.Ek_hist = self.model.Ek
        self.Ep_hist = self.model.Ep
        self.Wp_hist = self.model.Wp
        self.Em_hist = self.model.Em
        self.Et_hist = np.array(self.Em_hist) + np.array(self.Wp_hist)

        if self.animate:
            self.setup_animation()
    
    def spring_shape(self, y_mass, x, n_coils=12, amplitude=0.3):
   
        y_top = 0
        ys = np.linspace(y_top, y_mass, n_coils*2)
        xs = []

        for i in range(len(ys)):
            if i == 0 or i == len(ys)-1:
                xs.append(x)
            else:
                xs.append(amplitude * (-1)**i + x)

        return xs, ys

    def setup_animation(self):

        self.fig, self.ax  = plt.subplots(1, 3, figsize = (14, 6), tight_layout = True)

        self.line, = self.ax[0].plot([],[], 'k-', lw = 2)
        self.point, = self.ax[0].plot([],[], 'ro', markersize = 10)
        
        self.ax[0].axhline(0, color='black', lw=0.5)
        self.ax[0].set_xlabel("x (m)")
        self.ax[0].set_ylabel("y (m)")

        self.line1, = self.ax[1].plot([],[], 'b-', lw = 2)

        self.ax[1].axhline(0, color = 'black', lw = 0.5)
        self.ax[1].set_ylim(- max(np.abs(self.q_hist) * 1.5), max(np.abs(self.q_hist) * 1.5))
        self.ax[1].set_xlim(0, self.t_max)
        self.ax[1].set_xlabel('t (s)')

        if self.system == 'pendulum':

            self.ax[0].set_xlim(-self.L*1.2, self.L*1.2)
            self.ax[0].set_ylim(-self.L*1.2, self.L*1.2)

            self.ax[1].set_ylabel('x (m)')


        else: 

            self.ax[0].set_xlim(-1,1)
            self.ax[0].set_ylim(-1.2* max(np.abs(self.q_hist)), 1.2*max(np.abs(self.q_hist)))

            self.ax[1].set_ylabel('y (m)')
        

        self.lineEk, = self.ax[2].plot([],[], 'b-', lw=2, label="Ek")
        self.lineEp, = self.ax[2].plot([],[], 'r-', lw=2, label="Ep")
        self.lineWp, = self.ax[2].plot([],[], 'g-', lw=2, label="Wd")
        self.lineEm, = self.ax[2].plot([],[], 'k-', lw=2, label="Em")
        self.lineEt, = self.ax[2].plot([],[], 'y-', lw =2, label = 'Et')

        self.ax[2].set_xlim(0, self.t_max)
        self.ax[2].set_ylim(- np.abs(self.Em_hist[0]) * 1.5, np.abs(self.Em_hist[0]) * 1.5)
        self.ax[2].set_xlabel('t (s)')
        self.ax[2].set_ylabel('E (J)')
        self.ax[2].legend()
    
        #Aimation creation
        self.frame_step = 5
        self.anim = FuncAnimation(self.fig, self.update, frames = np.arange(0, len(self.t_hist), self.frame_step), interval = 10, blit = False, repeat = False)
        #write = PillowWriter(fps = 30)

        #if self.model == 'pendulum':
        #    self.anim.save('C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Damping_vibration\\figures\\pendulum.gif', write)
        #else:
        #    self.anim.save('C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Damping_vibration\\figures\\spring.gif', write)
        plt.show()

    def update(self, i):
        #Update lines

        if self.system == 'pendulum':

            x = self.L * np.sin(self.q_hist[i])
            y = - self.L * np.cos(self.q_hist[i])

            self.point.set_data([x], [y])
            self.line.set_data([0, x], [0 ,y])

            self.line1.set_data(self.t_hist[:i+1], self.q_hist[:i+1])
        else:
            x = 0 
            xs, ys = self.spring_shape(self.q_hist[i], x)

            self.point.set_data([x], [self.q_hist[i]])
            self.line.set_data(xs,  ys)

            self.line1.set_data(self.t_hist[:i+1], self.q_hist[:i+1])
            


        self.lineEk.set_data(self.t_hist[:i+1], self.Ek_hist[:i+1])
        self.lineEp.set_data(self.t_hist[:i+1], self.Ep_hist[:i+1])
        self.lineWp.set_data(self.t_hist[:i+1], self.Wp_hist[:i+1])
        self.lineEm.set_data(self.t_hist[:i+1], self.Em_hist[:i+1])
        self.lineEt.set_data(self.t_hist[:i+1], self.Et_hist[:i+1])

class Pendulum():

    g = 9.81

    def __init__(self, theta0, omega0, m, L, gamma, dt, t_max, approx=False):

        self.theta = theta0
        self.omega = omega0

        self.m = m
        self.L = L
        self.gamma = gamma

        self.dt = dt
        self.t_max = t_max


        if approx == True and np.rad2deg(theta0) > 15:
            raise ValueError('The initial angle must be less than 15')
        else:
            self.approx = approx

        # damping
        self.beta = gamma / (m)

        # histories
        self.t = 0
        self.t_hist = [0]
        self.theta_hist = [theta0]
        self.omega_hist = [omega0]

        self.Ek = []
        self.Ep = []
        self.Wp = []
        self.Em = []

        Ek, Ep, Wp, Em = self.energy(theta0, omega0, dt)

        self.Ek.append(Ek)
        self.Ep.append(Ep)
        self.Wp.append(Wp)
        self.Em.append(Em)
       


    def solve(self):

        while self.t < self.t_max:

            if self.approx:
                theta_new, omega_new = self.euler()

            elif self.beta > 5:
                theta_new, omega_new = self.crank_nicolson(self.dt)

            else:
                theta_new, omega_new = self.rk4()

            Ek, Ep, Wp, Em = self.energy(self.theta, self.omega, self.dt)

            self.Ek.append(Ek)
            self.Ep.append(Ep)
            self.Wp.append(Wp + self.Wp[-1])
            self.Em.append(Em)

            self.theta = theta_new
            self.omega = omega_new

            self.t += self.dt

            self.t_hist.append(self.t)
            self.theta_hist.append(self.theta)
            self.omega_hist.append(self.omega)

        self.q_hist = self.theta_hist
        self.dq_hist = self.omega_hist

    def euler(self):

        if abs(self.theta) > np.deg2rad(15):
            raise ValueError("Small-angle approximation violated")

        alpha = (
            -self.beta * self.omega
            - (self.g / self.L) * self.theta
        )

        omega_new = self.omega + alpha * self.dt
        theta_new = self.theta + self.omega * self.dt

        return theta_new, omega_new
    
    def rk4(self):

        dt = self.dt

        def f_theta(theta, omega):
            return omega

        def f_omega(theta, omega):
            return (-self.beta * omega - (self.g / self.L) * np.sin(theta))

        th = self.theta
        om = self.omega

        k1_th = f_theta(th, om)
        k1_om = f_omega(th, om)

        k2_th = f_theta(th + 0.5*dt*k1_th, om + 0.5*dt*k1_om)
        k2_om = f_omega(th + 0.5*dt*k1_th, om + 0.5*dt*k1_om)

        k3_th = f_theta(th + 0.5*dt*k2_th, om + 0.5*dt*k2_om)
        k3_om = f_omega(th + 0.5*dt*k2_th, om + 0.5*dt*k2_om)

        k4_th = f_theta(th + dt*k3_th, om + dt*k3_om)
        k4_om = f_omega(th + dt*k3_th, om + dt*k3_om)

        theta_new = th + dt/6*(k1_th+2*k2_th+2*k3_th+k4_th)
        omega_new = om + dt/6*(k1_om+2*k2_om+2*k3_om+k4_om)

        return theta_new, omega_new
    
    def crank_nicolson(self, dt):

        x_new = self.theta
        v_new = self.omega

        max_iter = 10
        tol = 1e-6

        for k in range(max_iter):

            F1 = x_new - self.theta - dt/2*(self.omega + v_new) 
            F2 = v_new - self.omega - dt/2*(-self.beta * (self.omega + v_new) - (self.g/self.L) * (np.sin(self.theta) + np.sin(x_new))) 

            # Jacobian
            dF1_dx = 1 
            dF1_dv = -dt/2 
            dF2_dx = -dt/2 * (-(self.g/self.L) * np.cos(x_new)) 
            dF2_dv = 1 - dt/2 * (-self.beta)
            
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
    
    def  energy(self, theta, omega, dt):
        # Disipacion de la energia por el teorema de las fuerzas vivas
        Ek = 0.5 * self.m * (omega*self.L) **2
        Ep = self.m * self.g * self.L *(1 - np.cos(theta))
        Wp = self.gamma * (omega* self.L) **2 * dt
        Em = Ek + Ep

        return Ek, Ep, Wp, Em
            
class Spring():

    g = 9.81

    def __init__(self, y0, v0, m, k, gamma, dt, t_max):

        # state variables
        self.y = y0
        self.v = v0

        self.m = m
        self.k = k
        self.gamma = gamma

        self.dt = dt
        self.t_max = t_max
        self.t = 0

        # parameters
        self.beta = gamma / (2*m)
        self.omega02 = k / m #Square

        # histories
        self.t_hist = [0]
        self.y_hist = [y0]
        self.v_hist = [v0]

        self.Ek = []
        self.Ep = []
        self.Wp = []
        self.Em = []

        Ek, Ep, Wp, Em = self.energy(y0, v0, self.dt)

        self.Ek.append(Ek)
        self.Ep.append(Ep)
        self.Wp.append(Wp)
        self.Em.append(Em)

    def solve(self):

        while self.t_max > self.t:

            y_new, v_new = self.rk4()
            
            Ek, Ep, Wp, Em = self.energy(y_new,  v_new, self.dt)

            self.Ek.append(Ek)
            self.Ep.append(Ep)
            self.Wp.append(Wp + self.Wp[-1])
            self.Em.append(Em)

            # Asign current state
            self.y = y_new
            self.v = v_new

            # Update current postion
            self.t += self.dt
            self.t_hist.append(self.t)
            self.y_hist.append(y_new)
            self.v_hist.append(v_new)

        self.q_hist = self.y_hist
        self.dq_hist = self.v_hist

    def rk4(self):

        dt = self.dt

        def f_y(y, v):
            return v

        def f_v(y, v):
            return -2*self.beta*v - self.omega02*y

        y = self.y
        v = self.v

        k1_y = f_y(y, v)
        k1_v = f_v(y, v)

        k2_y = f_y(y + 0.5*dt*k1_y, v + 0.5*dt*k1_v)
        k2_v = f_v(y + 0.5*dt*k1_y, v + 0.5*dt*k1_v)

        k3_y = f_y(y + 0.5*dt*k2_y, v + 0.5*dt*k2_v)
        k3_v = f_v(y + 0.5*dt*k2_y, v + 0.5*dt*k2_v)

        k4_y = f_y(y + dt*k3_y, v + dt*k3_v)
        k4_v = f_v(y + dt*k3_y, v + dt*k3_v)

        y_new = y + dt/6* (k1_y+2*k2_y+2*k3_y+k4_y)
        v_new = v + dt/6 *(k1_v+2*k2_v+2*k3_v+k4_v)

        return y_new, v_new
    
    def energy(self, y, v, dt):

        Ek = 0.5 * self.m * v**2
        Ep = 0.5 * self.k * y **2
        Wp = self.gamma * v**2 * dt
        Em = Ek + Ep 
        return Ek, Ep, Wp, Em


if __name__ == '__main__':
    sim = DampedVibration(q0=-2, dq0=0, m=1, gamma=0.3, t_max=10, dt=0.01, system="spring", animate=True, k=1)
    sim.run()