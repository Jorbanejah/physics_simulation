import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
"""
This code store in different .npz file the following data: 
- Lypunov coefficient through different values. (data_bifurcation.npz)
- Bifurcation diagram: alpha vs theta, alpaha vs dtheta. (data_bifurcation.npz)
- The trajectories of 1, 2, 4 - periods and chaotic motion. (data_poin_traj.npz)
- Poincare sections of chaotic motion. (poincare_sections.npz)

The following functions are commented and described. 
"""
@dataclass
class Params:
    beta: float = 0.5         # damping
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
    """
    Calculus of bifurcation diagram data

    Ags:
        f: function
        p: function params
        alphas: varible to change
        T: period
        t_transient (default): 100
        t_steady (default): 50
    """
    results_A = []
    results_poin_theta = []
    results_poin_omega = []

    for i, A in enumerate(alphas):
        p.A = A
        progress = (i + 1) / len(alphas)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)

        print(rf"Bifurcation diagram: [{bar}]  {progress*100:5.1f}%   A = {A:.4f}", end="\r", flush=True)

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

    
    return results_A, results_poin_theta, results_poin_omega

def trajectory(f, p, A, T, n_periods=80, dt=0.01):

    p.A = A

    t_final = n_periods * T

    sol = solve_ivp(f, [0, t_final], p.y0, args=(p,), max_step=dt, dense_output=True)

    t_vals = np.linspace(0, t_final, int(t_final/dt))
    theta_vals, omega_vals = sol.sol(t_vals)

    theta_vals = (theta_vals + np.pi) % (2*np.pi) - np.pi

    t_wrapped =(t_vals % T)
    return theta_vals, omega_vals, t_wrapped

def poincare_sections(f, A, p, T, n_trans=400, n_points=300):

    p.A = A
    theta_list = []
    omega_list = []

    # 1) Transient
    sol = solve_ivp(f, [0, n_trans*T], p.y0, args=(p,), max_step=0.05, dense_output=True)
    y = sol.y[:, -1]
    t0 = sol.t[-1]

    # 2) Period-by-period integration with absolute time
    for _ in range(n_points):
        sol = solve_ivp(f, [t0, t0 + T], y, args=(p,), max_step=0.05, dense_output=True)
        y = sol.y[:, -1]
        t0 = sol.t[-1]

        theta = sol.sol(t0)[0]
        omega = sol.sol(t0)[1]

        theta = (theta + np.pi) % (2*np.pi) - np.pi

        theta_list.append(theta)
        omega_list.append(omega)

    return theta_list, omega_list, None


def lyapunov_exponent(f, p, A, steps=8000, dt=0.01, delta0=1e-8):

    p.A = A

    # 1) Attractor first
    T = 2*np.pi / p.omega_drive
    sol = solve_ivp(f, [0, 200*T], p.y0, args=(p,), max_step=0.01)
    x_attractor = sol.y[:, -1]

    # 2) Use that as initial condition
    x1 = x_attractor.copy()

    np.random.seed(42)
    v = np.random.normal(size=2)
    v /= np.linalg.norm(v)
    x2 = x1 + delta0 * v

    S = 0.0
    t = 0.0

    for _ in range(steps):
        sol1 = solve_ivp(f, [t, t+dt], x1, args=(p,), max_step=dt)
        sol2 = solve_ivp(f, [t, t+dt], x2, args=(p,), max_step=dt)

        x1 = sol1.y[:, -1]
        x2 = sol2.y[:, -1]

        diff = x2 - x1
        dist = np.linalg.norm(diff)

        S += np.log(dist / delta0)

        diff = diff * (delta0 / dist)
        x2 = x1 + diff

        t += dt

    return S / (steps * dt)


def compute_and_store(alphas, filename="C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz"):

    p = Params()

    # 1) Bifurcation diagram + Lypaunov coefficient -------> ~ 5/7 min

    data_L_bi = {
        "alphas": None,
        "bifur_q": None,
        "bifur_dq": None,
        "Lyapunov": []
    }
    T = 2 * np.pi /p.omega_drive
       
    #results_A, q_b, dq_b= bifurcation_diagram(driven_equation, p, alphas, T, t_transient=100, t_steady=50)

    #data_L_bi["alphas"] = results_A
    #data_L_bi["bifur_q"] = q_b
    #data_L_bi["bifur_dq"] = dq_b
    
    #for i, A in enumerate(alphas):

    #    progress = (i + 1) / len(alphas)
    #    bar_length = 12
    #    filled = int(progress * bar_length)
    #    bar = "█" * filled + "-" * (bar_length - filled)

    #    print(rf"Lyapunov: [{bar}]  {progress*100:5.1f}%   A = {A:.4f}", end="\r", flush=True)

    #    Lam = lyapunov_exponent(driven_equation, p, A)

    #    data_L_bi["Lyapunov"].append(Lam)

    #data_L_bi["Lyapunov"] = np.array(data_L_bi["Lyapunov"])

    # Stored Bifurcation diagram and Lyapunov coefficient
    #np.savez(filename, **data_L_bi)

    # 2) Poincare sections (chaotic regime)
    
    data_P = {
        "poincare_q": {},
        "poincare_dq": {}
    }
    initial_y0 = np.array([0.0, 0.0])

    for i, alpha in enumerate(alphas):

        progress = (i + 1) / len(alphas)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)

        print(rf"Poincare sections: [{bar}]  {progress*100:5.1f}% $\alpha= {alpha:.4f}$", end="\r", flush=True)

        p.y0 = initial_y0.copy()
        theta_list, omega_list, _ = poincare_sections(driven_equation, alpha, p, T)

        data_P["poincare_q"][alpha] = theta_list
        data_P["poincare_dq"][alpha] = omega_list
    
    np.savez(filename, **data_P)

    # 3) Trajectories  + poincare sections ----> (~ 1 min)

    data_t_p = {
        "q_poincare": {},
        "dq_poincare": {},
        "q_trajectory": {},
        "dq_trajectory":{},
        "t_mod_T": {}
    }

    
    #for i, alpha in enumerate(A_values):

    #    progress = (i + 1) / len(alphas)
    #    bar_length = 12
    #    filled = int(progress * bar_length)
    #    bar = "█" * filled + "-" * (bar_length - filled)

    #    print(f"Poincare sections: [{bar}]  {progress*100:5.1f}% $\alpha$ = {alpha:.4f}", end="\r", flush=True)

    #    theta_vals, omega_vals, t_wrapped = trajectory(driven_equation, p, alpha, T)

    #    data_t_p["q_trajectory"][alpha] = theta_vals
    #    data_t_p["dq_trajectory"][alpha] = omega_vals
    #    data_t_p["t_mod_T"][alpha] = t_wrapped

    #    theta_list, omega_list, _ = poincare_sections(driven_equation, alpha, p, T)

    #    data_t_p["q_poincare"][alpha] = theta_list
    #    data_t_p["dq_poincare"][alpha] = omega_list

    #np.savez(filename, **data_t_p)
    
    print(f"\nSaved to {filename}")

    


from mpl_toolkits.mplot3d import Axes3D

def plot_trajectory_and_poincare(f, p, A_values):

    T = 2 * np.pi / p.omega_drive

    fig = plt.figure(figsize=(14, 12))
    
    # Bright, contrasting colors
    traj_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd"]
    poin_colors = ["#d62728", "#17becf", "#e377c2", "#8c564b"]

    for idx, A in enumerate(A_values):

        progress = (idx + 1) / len(A_values)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)
        print(rf"Poincare and trajectories: [{bar}]  {progress*100:5.1f}%   A = {A:.4f}",
              end="\r", flush=True)

        # 3D subplot
        ax = fig.add_subplot(2, 2, idx+1, projection='3d')

        # Reset initial condition
        p.y0 = np.array([0.0, 0.0])

        # Compute 3D trajectory
        theta_traj, omega_traj, t_wrap = trajectory(f, p, A, T, n_periods=300, dt=0.01)

        # Compute Poincaré section
        theta_poin, omega_poin, t_poin = poincare_sections(f, A, p, T)

        # Plot trajectory (3D curve)
        ax.plot(theta_traj, omega_traj, t_wrap,
                color=traj_colors[idx],
                alpha=0.55,
                linewidth=1.2,
                label="Trajectory")

        # Plot Poincaré points (all at same t mod T)
        ax.scatter(theta_poin, omega_poin, t_poin,
                   color=poin_colors[idx],
                   s=35,
                   edgecolor="black",
                   linewidth=0.4,
                   label="Poincaré")

        ax.set_title(f"A = {A}", fontsize=14, weight="bold")
        ax.set_xlabel(r"$\theta$")
        ax.set_ylabel(r"$\dot{\theta}$")
        ax.set_zlabel(r"$t mod(T)$")
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.suptitle("Driven Pendulum: 3D Trajectories and Poincaré Sections", fontsize=18, weight="bold")
    plt.show()

if __name__ == "__main__":
    p =Params()
    #alphas = np.linspace(1.060, 1.087, 100) #----> Bifurcation and Lyapunov coefficient 
    
    #A_values = [0.5, 1.07, 1.09, 1.5]  #-----> trajectories and poincare sections

    alphas = np.linspace(1.1, 1.6, 50) # -----> poincare sections (chaos regime)

    compute_and_store(alphas)

    




