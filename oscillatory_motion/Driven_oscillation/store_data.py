import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

# ============================================================
# 1. PARAMETERS AND EQUATION (Taylors form)
# ============================================================

@dataclass
class Params:
    beta: float = 0.5          # damping
    omega: float = 2/3          # drive frequency
    dt: float = 0.05            # time step
    A: float = 1.083
    y0: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0]))  # Fixed dataclass default



def driven_equation(t, y, p):
    theta, omega = y
    return np.array([omega, -p.beta * omega - np.sin(theta) + p.A * np.cos(p.omega* t)])

#=============================================================
# 2. INTEGRATOR
#=============================================================

def rk4(f, t, y, dt, params=None):
    """
    Generic RK4 integrator

    Parameters
    ----------
    f : function
        f(t, y, params) -> dy/dt
    t : float
    y : np.ndarray
    dt : float
    params : optional

    Returns
    -------
    y_new : np.ndarray
    """
    if params is None:
        raise ValueError("Params must be provided")

    k1 = f(t, y, params)
    k2 = f(t + dt/2, y + dt/2 * k1, params)
    k3 = f(t + dt/2, y + dt/2 * k2, params)
    k4 = f(t + dt, y + dt * k3, params)

    return y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

# ============================================================
# 3. POINCARÉ SECTION
# ============================================================

def poincare_section(p, A):
    """
    Compute Poincaré section by sampling once per drive period
    after a transient. Phase-locked to standard driven pendulum 
    bifurcation diagrams.
    
    Parameters
    ----------
    p : Params
        Parameter dataclass
    A : float
        Driving amplitude
    
    Returns
    -------
    theta_wrapped : np.ndarray
        Wrapped angles at Poincaré section
    omega_samples : np.ndarray  
        Angular velocities at Poincaré section
    results_A : np.ndarray
        Corresponding A values
    """
    #Calculus period
    T = 2.0 * np.pi / p.omega 

    #Define the total cycles
    t_transient = 100
    t_steady = 50

    t_total = (t_transient + t_steady) * T
    dt = p.dt

    y0 = p.y0.copy()
    t = 0.0

    p.A = A

    # 1) Compute full trajectory
    trajectory = []

    while t <= t_total:
        # 1) Total trajectory
        y = rk4(driven_equation, t = t, y = y0, dt = dt, params= p)
        trajectory.append(y)
        t += dt

    trajectory = np.array(trajectory)

    # 2) Poincaré sampling at exact multiples of T after transient

    t_eval = np.arange(t_transient * T, t_total, T) / dt  # Time indices
    t_eval = t_eval.astype(int)  # Integer indices for array slicing
    
    y_samples = trajectory[t_eval]
    theta_samples = y_samples[:, 0]
    omega_samples = y_samples[:, 1]
    
    # Wrap angle to [-π, π]
    theta_wrapped = (theta_samples + np.pi) % (2 * np.pi) - np.pi
    
    # Store results
    results_A = np.full(len(theta_wrapped), A)
    
    # Update initial condition for next A value (continuity)
    p.y0 = y_samples[-1].copy() 
    
    return theta_wrapped, omega_samples, results_A

# ============================================================
# 4. LYAPUNOV EXPONENT
# ============================================================

def lyapunov(p, x0, delta0=1e-8, h=0.01, t_final=100.0, renorm_interval=0.1): 
    
    x1 = np.array(x0, dtype=float) #Initial condition
    
    v = np.random.normal(size=3) #Random velocity 
    v /= np.linalg.norm(v) #Normalize
    x2 = x1 + delta0 * v #Initial condition perturbation

    S = 0.0  # Acumulative factor 
    t = 0.0
    steps_per_renorm = int(renorm_interval / h) 
    total_steps = int(t_final / h) #Total number of steps 
    
    #We are going to compute those A values that split the graphics (1 -T, 2- T, 3-T, chaotic)
    S_compute =[]
    As = 0
    full_trajectories = {A: {"q_full": {}, "d_full": {}} for A in As}

    for i in range(0, total_steps):
        x1 = rk4(driven_equation, t = i, dt = h, params = p)
        x2 = rk4(driven_equation, x2, h)
        if i % steps_per_renorm == 0: #El % en el condicional es el operador módulo, que devuelve el resto de la división
            diff = x2 - x1
            dist = np.linalg.norm(diff)
            S += np.log(dist / delta0)
            diff = (delta0 / dist) * diff
            x2 = x1 + diff
            t += renorm_interval
    S_compute.extend(S / (t))


# ============================================================
# 5. SWEEP OVER ALPHAS AND STORE EVERYTHING
# ============================================================

def compute_and_store(alphas, filename="C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\store_data.npz"):
    p = Params()
    data = {
        "q_poincare": {},
        "dq_poincare": {},
    }

    # 1) Poincare sections
    for i, A in enumerate(alphas):
       
        progress = (i + 1) / len(alphas)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)

        print(rf"[{bar}]  {progress*100:5.1f}%   $A$ = {A:.4f}", end="\r", flush=True)

        q_p, dq_p, _= poincare_section(p, A)

        # store
        data["q_poincare"][A] = q_p
        data["dq_poincare"][A] = dq_p

    np.savez(filename, **data)
    print(f"\nSaved to {filename}")


# ============================================================
# 6. RUN
# ============================================================

if __name__ == "__main__":

    alphas = np.linspace(1.060, 1.087, 50)
    compute_and_store(alphas, "store_data.npz")

