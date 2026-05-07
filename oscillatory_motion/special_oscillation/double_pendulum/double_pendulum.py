'''
The following code describe the motion of the double-pendulum system.
Description of the system: two masses join to a inextensable rod. The system will be need: L_1, L_2 -- string: m1, m2 -- masses
Structure of the code:

    Class Params: where we define the main parameters, on top of that a boolean parameter called: small-angle.

    def equation(params)
    def approx_equation(params)
    Class DoublePendulum: 
        def __init__(self, params, **kwargs):
        def run
        def Transform
        def energies
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass

@dataclass
class Params:
    g: float = 9.81  # m/s^2
    
    # Main parameters
    m1: float = 1.0  # kg
    m2: float = 1.5  # kg
    L1: float = 1.0  # m
    L2: float = 2.0  # m
    
    # Initial conditions
    q0: tuple = (np.deg2rad(30), 0)  # rad
    dq0: tuple = (0, 0)  # rad/s
    
    # Times
    t: float = 15.0  # s
    dt: float = 0.01

def velocity_verlet(f, t, y, dt, params=None):
    """"Velocity Verlet integrator for second-order system"""
    th1, w1, th2, w2 = y
    
    _, _, a1, a2 = f(t, y, p = params)
    
    th1_new = th1 + w1*dt + 0.5*a1*dt*dt
    th2_new = th2 + w2*dt + 0.5*a2*dt*dt

    _, _, a1_new, a2_new = f(t + dt, [th1_new, th2_new, w1, w2], params)

    w1_new = w1 + 0.5*(a1 + a1_new)*dt
    w2_new = w2 + 0.5*(a2 + a2_new)*dt

    return np.array([th1_new, th2_new, w1_new, w2_new])

def rk4(f, t, dt, y, p):
    """
    Generic RK4 integrator
    """
    k1 = f(t, y, p)
    k2 = f(t + dt/2, y + dt/2 * k1, p)
    k3 = f(t + dt/2, y + dt/2 * k2, p)
    k4 = f(t + dt, y + dt * k3, p)

    y_new = y + dt/6.0 * (k1 + 2 *k2 + 2* k3+ k4)

    return y_new

def equation_double_pendulum(t, y, p):
    """
    Full nonlinear double pendulum equations
    y = [th1, th2, w1, w2]
    """
    th1, th2, w1, w2 = y
    delta = th1 - th2
    
    m1, m2 = p.m1, p.m2
    L1, L2 = p.L1, p.L2
    g = p.g
    
    den_1 = L1 * (m1 + m2 * np.sin(delta)**2)
    den_2 = L2/L1 * den_1
    
    dw1 = (m2 * L1 * w1**2 * np.sin(delta) * np.cos(delta) +
           m2 * g * np.sin(th2) * np.cos(delta) +
           m2 * L2 * w2**2 * np.sin(delta) -
           (m1 + m2) * g * np.sin(th1)) / den_1
    
    dw2 = (-m2 * L2 * w2**2 * np.sin(delta) * np.cos(delta) +
           (m1 + m2) * g * np.sin(th1) * np.cos(delta) -
           (m1 + m2) * L1 * w1**2 * np.sin(delta) -
           (m1 + m2) * g * np.sin(th2)) / den_2
    
    return np.array([w1, w2, dw1, dw2])


def equation_approx(t, y, p):
    """
    Small angle approximation (linearized)
    sin(theta) = theta, cos(theta) = 1
    y = [th1, th2, w1, w2]
    """
    th1, th2, w1, w2 = y

    delta = th1 - th2
    m1, m2 = p.m1, p.m2
    L1, L2 = p.L1, p.L2
    g = p.g
    
    # Linearized equations for small angles
    denom1 = L1 * (2*m1 + m2 - m2*np.cos(delta)**2)

    d2th1 = (-g*(2*m1 + m2)*th1 + m2*g*np.sin(delta) + 
             2*m2*L2*w2**2*np.sin(delta)*np.cos(delta) + 
             2*m1*L1*w1**2*np.sin(delta)*np.cos(delta)) / denom1
    
    denom2 = L2 * (2*m1 + m2 - m2*np.cos(delta)**2)
    d2th2 = (2*(m1 + m2)*L1*w1**2*np.sin(delta) + 
             g*(m1 + m2)*(th1 - th2) + 
             m2*L2*w2**2*np.sin(delta)*np.cos(delta)) / denom2
    
    
    return np.array([w1, w2, d2th1, d2th2])


class DoublePendulum:
    """
    Double Pendulum system.

    Parameters
    ------
    m1 : float
        first mass pendulum.
    m2 : float
        second mass pendulum.
    L1 : float
        length first rod pendulum. 
    L2 : float
        length second rod pendulum.
    q0 : float
        Initial generalize coordenate (th1, th2)
    dq0 : float
        Initial generalize velocity (w1, w2)
    t: int, optional
        total system time (default: 15)
    dt: float, optional
        Time step for numerical integration (default: 0.01)

    Other parameters:
    ------

    Small_angle : bool, optional
        parameter that decides which equation use either non-approximation or approximation (default: False)

    Note:
    ------
    This class implements a double pendulum model that can operate in 
    either small-angle approximation or without approximation, depending on
    the selected configuration of small_angle paramter

    """
    
    def __init__(self, params: Params, small_angle: bool = False, method: str = 'RK4'):
        self.params = params
        self.small_angle = small_angle
        self.method = method
        self.sol = None
        self.t = None
        self.y = None
        
    def run(self):
        """Run simulation using specified method"""
        y0 = np.array([self.params.q0[0], self.params.q0[1], 
                       self.params.dq0[0], self.params.dq0[1]])
        
        if self.method == 'RK4':
            # RK4 implementation
            n_steps = int(self.params.t / self.params.dt)
            t_points = np.linspace(0, self.params.t, n_steps + 1)
            
            y_history = np.zeros((4, n_steps + 1))
            t_history = np.zeros(n_steps + 1)
            
            y_history[:, 0] = y0
            t_history[0] = 0
            
            current_y = y0.copy()
            current_t = 0
            
            for i in range(n_steps):
                if self.small_angle:
                    current_y = rk4(equation_approx, current_t, current_y, 
                                   self.params.dt, self.params)
                else:
                    current_y = rk4(equation_double_pendulum, current_t, current_y, 
                                   self.params.dt, self.params)
                
                current_t += self.params.dt
                y_history[:, i + 1] = current_y
                t_history[i + 1] = current_t
            
            self.sol = {'t': t_history, 'y': y_history}
            self.t = t_history
            self.y = y_history

        elif self.method =="Verlet":

            n_steps = int(self.params.t / self.params.dt)
            t_history = np.linspace(0, self.params.t, n_steps + 1)
            
            y_history = np.zeros((4, n_steps + 1))
            y_history[:, 0] = y0
            
            current_y = y0.copy()
            current_t = 0
            
            eq_func = equation_approx if self.small_angle else equation_double_pendulum
            
            for i in range(n_steps):
                current_y = velocity_verlet(eq_func, current_t, current_y, 
                                          self.params.dt, self.params)
                current_t += self.params.dt
                y_history[:, i + 1] = current_y
            
            self.sol = {'t': t_history, 'y': y_history}
            self.t = t_history
            self.y = y_history

        else:  # Default to solve_ivp RK45
            t_span = (0, self.params.t)
            t_eval = np.arange(0, self.params.t + self.params.dt/2, self.params.dt)
            
            if self.small_angle:
                sol = solve_ivp(equation_approx, t_span, y0, args=(self.params,), 
                              t_eval=t_eval, method='RK45', rtol=1e-8)
            else:
                sol = solve_ivp(equation_double_pendulum, t_span, y0, args=(self.params,), 
                              t_eval=t_eval, method='RK45', rtol=1e-8)
            
            self.sol = sol
            self.t = sol.t
            self.y = sol.y
        
        return self.sol
    
    def transform(self):

        """Convert to Cartesian coordinates"""
        if self.y is None:
            raise ValueError("Run simulation first!")
            
        th1, th2, w1, w2 = self.y
        
        # Position of first mass
        x1 = self.params.L1 * np.sin(th1)
        y1 = -self.params.L1 * np.cos(th1)
        
        # Position of second mass
        x2 = x1 + self.params.L2 * np.sin(th2)
        y2 = y1 - self.params.L2 * np.cos(th2)
        
        return x1, y1, x2, y2
    
    def energies(self):
        """Calculate kinetic and potential energies"""
        if self.y is None:
            raise ValueError("Run simulation first!")
            
        th1, th2, w1, w2 = self.y
        
        m1, m2, L1, L2, g = self.params.m1, self.params.m2, self.params.L1, self.params.L2, self.params.g
        
        # Positions
        x1, y1, x2, y2 = self.transform()
        
        # Velocities
        vx1 = L1 * w1 * np.cos(th1)
        vy1 = L1 * w1 * np.sin(th1)
        vx2 = vx1 + L2 * w2 * np.cos(th2)
        vy2 = vy1 + L2 * w2 * np.sin(th2)
        
        # Kinetic energy
        T1 = 0.5 * m1 * (vx1**2 + vy1**2)
        T2 = 0.5 * m2 * (vx2**2 + vy2**2)
        T = T1 + T2
        
        # Potential energy
        V1 = m1 * g * y1
        V2 = m2 * g * y2
        V = V1 + V2
        
        return T, V, T + V


# Example usage and plotting
def plot_double_pendulum():
    params = Params()
    
    # Full nonlinear simulation
    pendulum_full = DoublePendulum(params, small_angle=False)
    pendulum_full.run()
    
    # Small angle approximation
    pendulum_approx = DoublePendulum(params, small_angle=True)
    pendulum_approx.run()
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Phase plot theta1
    axes[0,0].plot(pendulum_full.y[0], pendulum_full.y[2], label='Full')
    axes[0,0].plot(pendulum_approx.y[0], pendulum_approx.y[2], '--', label='Approx')
    axes[0,0].set_xlabel('θ₁ (rad)')
    axes[0,0].set_ylabel('ω₁ (rad/s)')
    axes[0,0].legend()
    axes[0,0].grid(True)
    
    # Theta1 vs time
    axes[0,1].plot(pendulum_full.t, np.rad2deg(pendulum_full.y[0]), label='Full')
    axes[0,1].plot(pendulum_approx.t, np.rad2deg(pendulum_approx.y[0]), '--', label='Approx')
    axes[0,1].set_xlabel('Time (s)')
    axes[0,1].set_ylabel('θ₁ (deg)')
    axes[0,1].legend()
    axes[0,1].grid(True)
    
    # Transform to Cartesian and plot trajectory
    x1, y1, x2, y2 = pendulum_full.transform()
    axes[1,0].plot(x1, y1, 'b-', alpha=0.7, label='Mass 1')
    axes[1,0].plot(x2, y2, 'r-', alpha=0.7, label='Mass 2')
    axes[1,0].plot(0, 0, 'ko')
    axes[1,0].set_xlabel('x (m)')
    axes[1,0].set_ylabel('y (m)')
    axes[1,0].set_title('Trajectory')
    axes[1,0].legend()
    axes[1,0].grid(True)
    axes[1,0].axis('equal')
    
    # Energies
    T, V, E = pendulum_full.energies()
    
    axes[1,1].plot(pendulum_full.t, T, label='Kinetic')
    axes[1,1].plot(pendulum_full.t, V, label='Potential')
    axes[1,1].plot(pendulum_full.t, E, 'k--', label='Total')
    axes[1,1].set_xlabel('Time (s)')
    axes[1,1].set_ylabel('Energy (J)')
    axes[1,1].legend()
    axes[1,1].grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_double_pendulum()

