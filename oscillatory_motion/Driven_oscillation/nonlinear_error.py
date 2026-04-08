import numpy as np
import matplotlib.pyplot as plt
from time import time
from Driven_oscillation import DrivenOscillation
from dataclasses import dataclass   

# ============================================================
#   Parameter container
# ============================================================
@dataclass
class DrivenOscillationParams(): 
        
        #Pendulum
        mass: float = 2
        gamma: float = 1
        L: float = 2

        # External forces
        omega: float = 2
        F0: float = 1
        F_external: str = 'cos'
        system: str = 'nonlinear'

        #Initial condition
        q0: float = np.deg2rad(30)
        dq0: float = 0 

        #Times
        dt: float = 0.01
        t_max: float = 20
        
# ============================================================
#   Simulation and time
# ============================================================

def run_simulation(params: DrivenOscillationParams):
    """
    Run a driven oscillation simulation and return (history, runtime).
    """
    start = time()

    osc = DrivenOscillation(q0 = params.q0, dq0 = params.dq0, m =params.mass, gamma=params.gamma, F0 = params.F0, omega= params.omega, t= params.t_max, dt = params.dt, system =params.system,L = params.L,F_external =params.F_external)

    model = osc.run()
    runtime = time() - start
    return model.history, runtime

# ============================================================
#   Energy + diagnostics from history
# ============================================================

def total_energy(history, method: str):
    """
    Compute total energy: mechanical + work.
    """
    Ek = np.array(history[method]["Ek"])
    Ep = np.array(history[method]["Ep"])
    Wd = np.array(history[method]["Wp_diss"])
    Wf = np.array(history[method]["Wp_drive"])
    return Ek + Ep + (Wd - Wf)


def energy_from_history(history, method: str):
    """
    Energy drift = max(E) - min(E) from a given history.
    """
    E = total_energy(history, method)
    return np.max(E),  np.min(E), E[0]


def max_amplitude_from_history(history, method: str):
    """
    Maximum absolute angle from a given history.
    """
    q = history[method]["q"]
    return np.max(np.abs(q))

# ============================================================
#   Plotting utilities
# ============================================================

def plot_energy_vs_initial_angle(theta, results_energy, methods):
    
    fig, ax = plt.subplots(1, len(methods), figsize = (12, 8), tight_layout = True)

    # Ensure ax is iterable even if only one method
    if len(methods) == 1:
        ax = [ax]

    for i, method in enumerate(methods):
        Emax = [results_energy[method]["E_max"][th] for th in theta]
        Emin = [results_energy[method]["E_min"][th] for th in theta]

        E0 = np.array([results_energy[method]["E0"][th] for th in theta])

        dE = (np.array(Emax) - np.array(Emin)) / E0

        ax[i].plot(theta, dE, label=r"$\Delta E / E_0$")

        ax[i].set_title(f"{method}")
        ax[i].set_xlabel("Initial angle (rad)")
        ax[i].set_ylabel(r"Energy drift $\Delta E / E_0$")
        ax[i].grid(alpha=0.3)
        ax[i].legend()

    fig.suptitle("Energy Drift vs Initial Angle", fontsize=14)
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\plot_energy_vs_initial_angle_{params.dt}.png", dpi=300, bbox_inches='tight')

def plot_runtime(dt_values, Errors):

    plt.figure(figsize=(6, 5))

    plt.plot(dt_values, Errors['rk4']["time"], "o-", lw=2)
    
    plt.xlabel(r"Time step $\Delta t$")
    plt.ylabel("Runtime (s)")
    plt.title("Runtime vs Time Step")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\runtime_nonlinear.png", dpi=300, bbox_inches='tight')

def plot_stability(dt_values, Errors, methods):

    fig, ax = plt.subplots(1, len(methods), figsize = (12, 8), tight_layout = True)

    if len(methods) == 1:
        ax = [ax]

    for i, method in enumerate(methods):
        ax[i].plot(dt_values, Errors[method]["max_amplitude"], "o-", lw=2)
        ax[i].set_title(f"Numerical Stability ({method})")
        ax[i].set_xlabel(r"Time step $\Delta t$")
        ax[i].set_ylabel(r"Max $|\theta|$")
        ax[i].grid(alpha=0.3)

    fig.suptitle("Stability vs Time Step", fontsize=14)
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\stability_nonlinear.png", dpi=300, bbox_inches='tight')

def plot_convergence(dt_values, Errors, methods):
    """
    Plot log(error) vs log(dt) and include theoretical slope lines.
    """
    dt_values = np.array(dt_values, dtype=float)

    fig, ax = plt.subplots(1, len(methods), figsize = (12, 8), tight_layout = True)

    if len(methods) == 1:
        ax = [ax]

    theoretical_order = {
            "CN": 2,
            "Verlet": 2
        }
    
    for i, method in enumerate(methods):
        err = np.array(Errors[method]["Error"])
        
        ax[i].loglog(dt_values, err, "o-", lw=2, markersize=6, label=method)

        # Add theoretical slope line if known
        order = theoretical_order.get(method, None)

        if order is not None:
            x0 = dt_values[len(dt_values) // 2]
            y0 = err[len(err) // 2]
            slope_line = y0 * (dt_values / x0) ** order
            ax[i].loglog(dt_values, slope_line, "--", color="gray", label=rf"Slope ≈ {order}")

        ax[i].set_title(f"Convergence ({method})")
        ax[i].set_xlabel(r"Time step $\Delta t$")
        ax[i].set_ylabel("Error")
        ax[i].grid(True, which="both", alpha=0.3)
        ax[i].legend()

    fig.suptitle("Convergence Test", fontsize=14)
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\convergence_nonlinear.png", dpi=300, bbox_inches='tight')

# ============================================================
#   Main error analysis
# ============================================================

def Energy_vs_initial_angle(params, methods):
    """
    ---------------------------------------------------------
    1. Energy drift vs initial angle
    ---------------------------------------------------------
    Computes the maximum and minimum energy for each numerical
    method as a function of the initial angle. For each angle,
    the simulation is run and the energy drift is extracted.
    Results are stored in a nested dictionary:

        results_energy[method]["E_max"][theta]
        results_energy[method]["E_min"][theta]

    Finally, the function plots energy drift vs initial angle.
    """
    print("Energy drift vs initial angle")

    results_energy = {
        name: {"E_max": {},
        "E_min": {}, "E0": {}} for name in methods
    }
    
    theta_deg = np.arange(1, 81)
    theta_values = np.deg2rad(theta_deg)
    params.dt = 0.001
    for th in theta_values:
        params.q0 = th
        history, _ = run_simulation(params)
        for method in methods:
            Emax, Emin, E0 = energy_from_history(history, method)
            results_energy[method]["E_max"][th] = Emax
            results_energy[method]["E_min"][th] = Emin
            results_energy[method]["E0"][th] = E0


    plot_energy_vs_initial_angle(theta_values, results_energy, methods)
    

def convergence_runtime(params, dt_values, methods):
    """
    ---------------------------------------------------------
    2. Convergence test (reference RK4) and runtime
    ---------------------------------------------------------
    For each timestep dt:
        - Runs the simulation
        - Uses RK4 as reference solution
        - Computes RMS error for each method
        - Stores runtime and max amplitude
    Produces:
        - Convergence plot
        - Runtime plot
        - Stability plot
    """
    
    Errors = {name: {"Error": [], "dt": [], "time":[], "max_amplitude": []} for name in methods}
    
    history, _ = run_simulation(params)
    #Reference RK4
  

    for dt in dt_values:
        params.dt = dt
        history, time = run_simulation(params)
        t_ref = history["rk4"]["t"]
        q_ref = history["rk4"]["q"]
        for method in methods:
            Errors[method]["dt"].append(dt)
            Errors[method]["time"].append(time)
            Errors[method]["max_amplitude"].append(max_amplitude_from_history(history, method))

            theta_num = history[method]["q"]
            t_num = history[method]["t"]

            #Interpolate reference onto numerical grid
            theta_ref_interp = np.interp(t_num, t_ref, q_ref)

            err = np.sqrt(np.mean((theta_num - theta_ref_interp) ** 2))
            Errors[method]["Error"].append(err)

    print("Plotting convergence")
    plot_convergence(dt_values, Errors, methods[1:])
    print("Plotting runtime")
    plot_runtime(dt_values, Errors)
    print("Plotting stability")
    plot_stability(dt_values, Errors, methods) # 

    plt.show()


if __name__ == "__main__":
    
    """
    Run all numerical analyses for the driven nonlinear oscillator.

    Produces:
        1. Energy drift vs initial angle
        2. Convergence test
        3. Numerical stability test
        4. Runtime vs dt
    """

    params = DrivenOscillationParams()
    dt_values = [0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
    methods = ['rk4', 'CN', 'Verlet']

    Energy_vs_initial_angle(params, methods)

    convergence_runtime(params, dt_values, methods)

    plt.show()