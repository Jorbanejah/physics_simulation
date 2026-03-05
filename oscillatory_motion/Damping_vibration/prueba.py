import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# ============================================================
# Integrators
# ============================================================

def rk4_step(f, state, dt):
    """
    Generic RK4 integrator for first-order ODE systems.
    """

    k1 = f(state)
    k2 = f(state + 0.5 * dt * k1)
    k3 = f(state + 0.5 * dt * k2)
    k4 = f(state + dt * k3)

    return state + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)


def euler_step(f, state, dt):
    """
    Simple Euler integrator.
    """

    return state + dt * f(state)


# ============================================================
# Base dynamical system
# ============================================================

class DynamicalSystem:
    """
    Base class for damped dynamical systems.
    """

    def __init__(self, q0, dq0, m, gamma):

        self.q = q0
        self.dq = dq0

        self.m = m
        self.gamma = gamma

        self.t = 0

        # histories
        self.t_hist = [0]
        self.q_hist = [q0]
        self.dq_hist = [dq0]

        self.Ek = []
        self.Ep = []
        self.Wp = []
        self.Em = []

    def derivatives(self, state):
        raise NotImplementedError

    def energy(self, q, dq, dt):
        raise NotImplementedError


# ============================================================
# Pendulum
# ============================================================

class Pendulum(DynamicalSystem):

    g = 9.81

    def __init__(self, theta0, omega0, m, L, gamma):

        super().__init__(theta0, omega0, m, gamma)

        self.L = L

        # initial energy
        Ek, Ep, Wp, Em = self.energy(theta0, omega0, 0)
        self.Ek.append(Ek)
        self.Ep.append(Ep)
        self.Wp.append(Wp)
        self.Em.append(Em)

    def derivatives(self, state):

        theta, omega = state

        dtheta = omega
        domega = -(self.gamma/self.m)*omega - (self.g/self.L)*np.sin(theta)

        return np.array([dtheta, domega])

    def energy(self, theta, omega, dt):

        Ek = 0.5 * self.m * (omega*self.L)**2
        Ep = self.m * self.g * self.L * (1 - np.cos(theta))

        Wp = - self.gamma * (omega*self.L)**2 * dt

        Em = Ek + Ep

        return Ek, Ep, Wp, Em


# ============================================================
# Spring oscillator
# ============================================================

class Spring(DynamicalSystem):

    def __init__(self, y0, v0, m, k, gamma):

        super().__init__(y0, v0, m, gamma)

        self.k = k

        Ek, Ep, Wp, Em = self.energy(y0, v0, 0)

        self.Ek.append(Ek)
        self.Ep.append(Ep)
        self.Wp.append(Wp)
        self.Em.append(Em)

    def derivatives(self, state):

        y, v = state

        dy = v
        dv = -(self.gamma/self.m)*v - (self.k/self.m)*y

        return np.array([dy, dv])

    def energy(self, y, v, dt):

        Ek = 0.5 * self.m * v**2
        Ep = 0.5 * self.k * y**2

        Wp = - self.gamma * v**2 * dt

        Em = Ek + Ep

        return Ek, Ep, Wp, Em


# ============================================================
# Simulation
# ============================================================

class Simulation:

    def __init__(self, system, dt=0.01, t_max=20, integrator="rk4"):

        self.system = system

        self.dt = dt
        self.t_max = t_max

        if integrator == "rk4":
            self.integrator = rk4_step
        else:
            self.integrator = euler_step

    def run(self):

        state = np.array([self.system.q, self.system.dq])

        while self.system.t < self.t_max:

            state = self.integrator(self.system.derivatives, state, self.dt)

            q, dq = state

            self.system.q = q
            self.system.dq = dq

            self.system.t += self.dt

            Ek, Ep, Wp, Em = self.system.energy(q, dq, self.dt)

            self.system.Ek.append(Ek)
            self.system.Ep.append(Ep)
            self.system.Wp.append(self.system.Wp[-1] + Wp)
            self.system.Em.append(Em)

            self.system.t_hist.append(self.system.t)
            self.system.q_hist.append(q)
            self.system.dq_hist.append(dq)

class Duffing(DynamicalSystem):

    def derivatives(self,state):

        x,v = state

        dx = v
        dv = -delta*v - alpha*x - beta*x**3

        return np.array([dx,dv])

# ============================================================
# Animation
# ============================================================

class Animator:

    def __init__(self, system):

        self.sys = system

        self.fig, self.ax = plt.subplots(1,3,figsize=(14,5))

        # mechanical system
        self.point, = self.ax[0].plot([],[],'ro',markersize=8)
        self.line, = self.ax[0].plot([],[],'k-')

        # coordinate evolution
        self.line_q, = self.ax[1].plot([],[],'b-')

        # energies
        self.lineEk, = self.ax[2].plot([],[],'b-',label="Ek")
        self.lineEp, = self.ax[2].plot([],[],'r-',label="Ep")
        self.lineEm, = self.ax[2].plot([],[],'k-',label="Em")

        self.ax[1].set_xlim(0,self.sys.t_hist[-1])
        self.ax[2].set_xlim(0,self.sys.t_hist[-1])

        self.ax[2].legend()

    def update(self,i):

        q = self.sys.q_hist[i]

        if isinstance(self.sys,Pendulum):

            x = self.sys.L*np.sin(q)
            y = -self.sys.L*np.cos(q)

            self.point.set_data([x],[y])
            self.line.set_data([0,x],[0,y])

        else:

            self.point.set_data([0],[q])
            self.line.set_data([0,0],[0,q])

        self.line_q.set_data(self.sys.t_hist[:i],self.sys.q_hist[:i])

        self.lineEk.set_data(self.sys.t_hist[:i],self.sys.Ek[:i])
        self.lineEp.set_data(self.sys.t_hist[:i],self.sys.Ep[:i])
        self.lineEm.set_data(self.sys.t_hist[:i],self.sys.Em[:i])

        return self.point,

    def animate(self):

        anim = FuncAnimation(
            self.fig,
            self.update,
            frames=len(self.sys.t_hist),
            interval=30
        )

        plt.show()


# ============================================================
# Example usage
# ============================================================

if __name__ == "__main__":

    # choose system

    pendulum = Pendulum(
        theta0=0.2,
        omega0=0,
        m=1,
        L=1,
        gamma=0.2
    )

    sim = Simulation(pendulum,dt=0.01,t_max=20)

    sim.run()

    Animator(pendulum).animate()