import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dataclasses import dataclass

# =====================================================
# INTEGRATORS
# =====================================================

def rk4_step(f, state, dt):
    k1 = f(state)
    k2 = f(state + 0.5*dt*k1)
    k3 = f(state + 0.5*dt*k2)
    k4 = f(state + dt*k3)

    return state + dt/6*(k1+2*k2+2*k3+k4)


def symplectic_euler(force, q, p, dt):
    """
    Semi-implicit Euler (symplectic)
    Stable for strong damping/stiff systems
    """
    p = p + dt * force(q, p)
    q = q + dt * p
    return q, p


# =====================================================
# BASE BODY
# =====================================================

class PhysicalSystem:

    def derivatives(self, state):
        raise NotImplementedError

    def energy(self, state):
        raise NotImplementedError

# =====================================================
# DAMPED PENDULUM
# =====================================================

@dataclass
class DampedPendulum(PhysicalSystem):

    m: float
    L: float
    gamma: float
    g: float = 9.81

    @property
    def beta(self):
        return self.gamma/(self.m)

    def derivatives(self, state):
        theta, omega = state

        d_theta = omega
        d_omega = -self.beta*omega - (self.g/self.L)*np.sin(theta)

        return np.array([d_theta, d_omega])

    def force(self, q, p):
        return -self.beta*p - (self.g/self.L)*np.sin(q)

    def energy(self, state):
        theta, omega = state

        Ek = 0.5*self.m*(self.L*omega)**2
        Ep = self.m * self.g * self.L * (1-np.cos(theta))

        return Ek, Ep, Ek+Ep


# =====================================================
# DAMPED SPRING
# =====================================================

@dataclass
class DampedSpring(PhysicalSystem):

    m: float
    k: float
    gamma: float

    @property
    def beta(self):
        return self.gamma/(2*self.m)

    def derivatives(self, state):
        x, v = state
        return np.array([v, -self.beta*v-(self.k/self.m)*x])

    def energy(self, state):
        x, v = state

        Ek = 0.5*self.m*v**2
        Ep = 0.5*self.k*x**2

        return Ek, Ep, Ek+Ep

class PhysicsEngine:

    def __init__(self, system, state0, dt):

        self.sys = system
        self.state = np.array(state0, float)
        self.dt = dt

        self.history = []
        self.energy = []

    def step(self):

        # AUTO INTEGRATOR SWITCH
        if hasattr(self.sys, "beta") and self.sys.beta > 5:
            q, p = self.state
            q, p = symplectic_euler(
                self.sys.force, q, p, self.dt
            )
            self.state = np.array([q, p])

        else:
            self.state = rk4_step(
                self.sys.derivatives,
                self.state,
                self.dt
            )

        self.history.append(self.state.copy())
        self.energy.append(
            self.sys.energy(self.state)
        )

def animate_system(engine, steps=2000):

    fig = plt.figure(figsize=(10,5))

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    line, = ax1.plot([], [], 'o-', lw=2)

    Ek_line, = ax2.plot([],[],label="Ek")
    Ep_line, = ax2.plot([],[],label="Ep")
    Em_line, = ax2.plot([],[],label="Em")

    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_aspect('equal')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    
    ax2.set_xlim(0, steps)
    ax2.set_ylim(0, engine.sys.energy(engine.state)[2] * 1.1)
    ax2.legend()

    tdata=[]
    Ek=[]
    Ep=[]
    Em=[]

    def update(frame):

        engine.step()

        theta, omega = engine.state

        x=np.sin(theta)
        y=-np.cos(theta)

        line.set_data([0,x],[0,y])

        ek,ep,em=engine.energy[-1]

        tdata.append(frame)
        Ek.append(ek)
        Ep.append(ep)
        Em.append(em)

        Ek_line.set_data(tdata,Ek)
        Ep_line.set_data(tdata,Ep)
        Em_line.set_data(tdata,Em)

        return line,

    ani = FuncAnimation(fig, update, frames=steps, interval=10)

    plt.show()

def phase_space_pendulum():

    thetas=np.linspace(0.1,np.pi,40)

    plt.figure()

    for θ0 in thetas:

        pend=DampedPendulum(1,1,0.2)
        eng=PhysicsEngine(pend,[θ0,0],0.01)

        for _ in range(800):
            eng.step()

        hist=np.array(eng.history)

        plt.plot(hist[:,0], hist[:,1], color=plt.cm.viridis(θ0/np.pi))

    plt.xlabel("Theta (rad) ")
    plt.ylabel("Omega (rad/s)")
    plt.title("Pendulum Phase Space")
    plt.show()

def phase_space_spring():

    ks=np.linspace(1,25,40)

    plt.figure()

    for k in ks:

        spring=DampedSpring(1,k,0.2)
        eng=PhysicsEngine(spring,[1,0],0.01)

        for _ in range(600):
            eng.step()

        hist=np.array(eng.history)

        plt.plot(hist[:,0],
                 hist[:,1],
                 color=plt.cm.plasma(k/25))

    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.title("Spring Phase Space")
    plt.show()


pend = DampedPendulum(m=1,L=1,gamma=8)

engine = PhysicsEngine(
    pend,
    state0=[np.pi/2,0],
    dt=0.01
)

animate_system(engine)

phase_space_pendulum()
phase_space_spring()