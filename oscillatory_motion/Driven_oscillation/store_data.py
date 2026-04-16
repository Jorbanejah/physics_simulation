import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field

@dataclass
class Params:
    beta: float = 0.5          # damping
    omega_drive: float = 2/3          # drive frequency
    A: float = 1.083
    y0: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0]))  

    def __post_init__(self):
        "Validation"
        assert self.beta > 0, "Damping must be positive"
        assert self.A > 0, "Drive amplitude must be positive"

def driven_equation(t, y, p):
    """
    Driven pendulum Taylor's form: theta'' + beta * theta' + sin(theta) = A cos(ω_drive * t)

    Args:
        t: time
        y: [theta, omega] state vector
        p: parameters
        
    Returns:
        [dtheta/dt, domega/dt]
    """
    theta, omega = y
    return [omega, -p.beta * omega - np.sin(theta) + p.A * np.cos(p.omega_drive * t)]

def bifurcation_diagram(f, p, alphas, T, t_transient = 100, t_steady = 50):

    results_A = []
    results_poin_theta = []
    results_poin_omega = []

    for i, A in enumerate(alphas):
        p.A = A
        progress = (i + 1) / len(alphas)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)

        print(rf"[{bar}]  {progress*100:5.1f}%   A = {A:.4f}", end="\r", flush=True)

        # Integrate full transient + steady window
        sol = solve_ivp(f, [0, (t_transient + t_steady) * T], p.y0, args=(p,), 
                        dense_output=True,   # IMPORTANT: allows exact sampling at multiples of T
                        max_step=0.05        # keeps integration stable
        )

        # Sample exactly at multiples of T
        t_eval = np.arange(t_transient * T, (t_transient + t_steady) * T, T)
        y_samples = sol.sol(t_eval)
        
        theta_samples = y_samples[0]
        omega_samples = y_samples[1]

        # Wrap angle to [-pi, pi]
        theta_wrapped = (theta_samples + np.pi) % (2*np.pi) - np.pi
        
        # Store results
        results_A.extend([A] * len(theta_wrapped))
        results_poin_theta.extend(theta_wrapped)
        results_poin_omega.extend(omega_samples)

        # Update initial condition for next A (follows attractor branch)
        p.y0 = sol.y[:, -1]

    
    return results_A, results_poin_omega, results_poin_theta

def poincare_sections(f, p, alphas, T, t_transient = 500, t_steady = 50):
    results_A = []
    results_poin_theta = []
    results_poin_omega = []

    for i, A in enumerate(alphas):
        p.A = A
        progress = (i + 1) / len(alphas)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)

        print(rf"[{bar}]  {progress*100:5.1f}%   A = {A:.4f}", end="\r", flush=True)

        # Integrate full transient + steady window
        sol = solve_ivp(f, [0, (t_transient + t_steady) * T], p.y0, args=(p,), 
                        dense_output=True,   # IMPORTANT: allows exact sampling at multiples of T
                        max_step=0.05        # keeps integration stable
        )

        # Sample exactly at multiples of T
        t_eval = np.arange(0, (t_transient + t_steady) * T, T)
        y_samples = sol.sol(t_eval)
        
        theta_samples = y_samples[0]
        omega_samples = y_samples[1]

        # Wrap angle to [-pi, pi]
        theta_wrapped = (theta_samples + np.pi) % (2*np.pi) - np.pi
        
        # Store results
        results_A.extend([A] * len(theta_wrapped))
        results_poin_theta.extend(theta_wrapped)
        results_poin_omega.extend(omega_samples)

        # Update initial condition for next A (follows attractor branch)
        p.y0 = sol.y[:, -1]

    
    return results_A, results_poin_omega, results_poin_theta
def compute_and_store(alphas, filename="C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\store_data.npz"):

    p = Params()

    data = {
        "full_trajectory": {},
        "results_alpha": [],
        "q_poincare": [],
        "dq_poincare": []
    }

    # 1) Poincare sections
    
    T = 2 * np.pi /p.omega
       
    results_A, dq_p, q_p, _= bifurcation_diagram(driven_equation, p, alphas, T, t_transient=100, t_steady=50)

    data["results_alpha"] = results_A
    data["q_poincare"] = q_p
    data["dq_poincare"] = dq_p

    np.savez(filename, **data)
    print(f"\nSaved to {filename}")
    
if __name__ == "__main__":
    p =Params()
    alphas = np.linspace(1.060, 1.087, 50)
    alpha = [1.087]
    
    A, dq, q = poincare_sections(driven_equation, p, alpha, T = 2*np.pi/p.omega_drive)

    #fig, (ax2 ) = plt.subplots(1, 2, figsize=(15, 8))
    x = np.sin(q)
    y = np.cos(q)
    plt.scatter(x,y, s = 10)
    
    plt.show()