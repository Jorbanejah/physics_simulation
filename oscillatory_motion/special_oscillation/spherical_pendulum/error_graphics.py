"""
Spherical pendulum - main graphics:

Numerical graphics can be run by both: linearized and normal simulation.

- Comparison method 
- Energy drift colormap 
- Runtime. 
- Convergence/stability.
- Density kde

"""

import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import pandas as pd

from spherical_pendulum import Spherical_Pendulum
from dataclasses import dataclass
from typing import Tuple, Sequence, Dict

def _as_pair(values: Sequence[float], name: str) -> Tuple[float, float]:
    if len(values) != 2:
        raise ValueError(f"{name} must contain exactly two values")
    return (float(values[0]), float(values[1]))

@dataclass
class Params:
    "Phisical parameters and initial condition for the pendulum"
    g: float = 9.81

    m: float = 1.0
    L:float = 2.0

    q0: tuple[float, float]= (np.deg2rad(45.0), np.deg2rad(10.0))
    dq0: tuple[float, float] = (0.0, 1.0)

    t:float = 15
    dt: float = 0.01

    def __post_init__(self) -> None:
        if self.g <= 0.0:
            raise ValueError("Gravity mus be positive.")
        if self.m <= 0.0:
            raise ValueError("Mass must be positive.")
        if self.L <= 0.0:
            raise ValueError("Lenghts must be positive.")
        
        self.q0 = _as_pair(self.q0, "q0")
        self.dq0 = _as_pair(self.dq0, "dq0")
           
##
# ------------------------- RUN INSTANCE -------------------------------
##

def run(dt: Sequence[float], time: int = 150, flag: bool = False,) -> Dict:
    """
    Compute energy drift and runtime for several integration methods.

    Parameters
    ----------
    dt : Sequence[float]
        Time-step values to test.
    time : int
        Maximum simulation time (default = 150).
    flag : bool
        Whether to use the linearized double pendulum (default= False).

    Returns
    -------
    Dict
        Nested dictionary: results[method_name][dt] = (result, analysis)
    """
    params = Params()

    methods_to_test = ["Rk4", "RK45", "Verlet", "DOP853"]

    #Large time
    params.t = time

    results_dict = {}

    #Run the all methods.
    for name in methods_to_test:

        if name == "Verlet" and flag == False:
            continue
        
        sim = Spherical_Pendulum(small_angle = flag, method = name)

        results_dict[name] = {}
        for time_step in dt:

            params.dt = time_step

            print(f"\nRunning with {name} with {time_step}...")
            result, runtime = sim.run(params = params)
            _, _, Et = sim.energies()
            drift_energy = np.abs(np.abs(max(Et)) - np.abs(min(Et)))

            results_dict[name][time_step] = (result, runtime, drift_energy)
    
    return results_dict

def upload_data(name:str = "stored_data.npz", directory: str = os.getcwd()) -> Dict:
    "Upload the grid variable"
    path = os.path.join(directory, name)

    if os.path.exists(path):
        
        data = np.load(path, allow_pickle= True)

        grid = data["grid"].item()

        return grid

    else:
        raise FileNotFoundError("Compute first the store_data.py file")
    

##
# ------------------------ NUMERICAL STUDY METHODS: time_compute, stability, convergence, kde_phase_space ---------
##

def time_compute(dt: Sequence[float], results_dict: Dict) -> plt.Figure:
    """
    Plot computation time vs time step.
    """
    fig = plt.figure(figsize=(10, 6))
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#ff7f0e']

    for i, (method_name, dt_dict) in enumerate(results_dict.items()):
        color = colors[i % len(colors)]
        times = []
        dts = []

        for time_step in sorted(dt_dict.keys()):
            
            if time_step not in dt:
                continue
            
            result, runtime, analysis = dt_dict[time_step]
            
            times.append(runtime)
            dts.append(time_step)

        # Sort for connected lines
        sorted_pairs = sorted(zip(dts, times))
        dts_sorted = [x for x, _ in sorted_pairs]
        times_sorted = [y for _, y in sorted_pairs]

        plt.plot(dts_sorted, times_sorted, 'o-', label=method_name, color = color, linewidth=2, markersize=6)
    
    plt.xlabel('Time Step ($dt$)')
    plt.ylabel('Run Time (s)')
    plt.title('Computational Cost vs Time Step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    
    return fig

def convergence(dt: Sequence[float], results_dict: Dict) -> plt.Figure:
    """
    Study numerical convergence: Error vs Time Step.
    """
    rows = []

    for method_name, dt_dict in results_dict.items():
        for time_step in dt:
            if time_step not in dt_dict:
                continue

            result, _, energy_drift = dt_dict[time_step]
            error = energy_drift 

            rows.append({
                "dt": time_step,
                "error": error,
                "method": method_name
            })

    df = pd.DataFrame(rows)

    fig = plt.figure(figsize=(10, 6))
    sns.lineplot(data=df, x="dt", y="error", hue="method", marker="o", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Time Step $dt$")
    plt.ylabel("Energy Drift (Error)")
    plt.ylim(0, 1e-4)
    plt.title("Convergence Study: |Error| vs Time Step")
    plt.grid(True, alpha=0.3, which="both")

    return fig

def stability(dt: Sequence[float], results_dict: Dict) -> plt.Figure:
    """
    Analyze numerical stability limits using phase space visualization.
    
    We check if the solution remains bounded or explodes (NaN/Inf).
    """
    rows = []
    
    for method_name, dt_dict in results_dict.items():
        for time_step in dt:

            if time_step not in dt_dict:
                continue
                
            result, _, energy_drift = dt_dict[time_step]
            
            # Check trajectory for NaNs or Infs
            E = energy_drift
           
            has_nan = np.any(np.isnan(E))
            has_inf = np.any(np.isinf(E))
            THRESHOLD = 1e2   # or 1e8, depends on your system
            has_blowup = np.any(np.abs(E) > THRESHOLD)
            
            # Determine stability label
            status = "Unstable" if (has_nan or has_inf or has_blowup) else "Stable"
            
            rows.append({
                "dt": time_step,
                "method": method_name,
                "status": status,
                "final_energy": E if not has_nan and not has_inf else np.nan
            })
            
    df = pd.DataFrame(rows)
    
    # Visualization: Scatter with stability coloring
    fig = plt.figure(figsize=(10, 6))
    
    # Map status to colors
    palette = {"Stable": "blue", "Unstable": "red"}
    
    sns.scatterplot(data=df, x="dt", y="method", hue="status", palette=palette, s=200, style="status")
    
    plt.title("Numerical Stability Map")
    plt.xlabel("Time Step $dt$")
    plt.ylabel("Integration Method")
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    
    return fig

def kde_phase_space_subplots(results_dict: Dict) -> plt.Figure:

    """
    Phase-space density KDE plots for a single integration method.
    Expects: results_dict[dt] = (result, analysis)
    """

    # ---- Choose smallest dt ----
    time_step = min(results_dict["DOP853"].keys())
    result, _, analysis = results_dict["DOP853"][time_step]

    # ---- Extract solution arrays ----
    # Case 1: result object with attributes
    if hasattr(result, "theta"):
        sol = {
            "theta": result.theta1,
            "dtheta": result.theta2,
            "phi": result.omega1,
            "dphi": result.omega2,
        }

    # Case 2: SciPy OdeResult (result.y is array)
    elif hasattr(result, "y"):
        idx_map = {"theta": 0, "dtheta": 1, "phi": 2, "dphi": 3}
        sol = {name: result.y[idx] for name, idx in idx_map.items()}

    else:
        raise ValueError("Result object has no recognizable structure.")

    # ---- Variables to plot ----
    var_pairs = [
        ("theta", "dtheta"),
        ("theta", "dphi"),
        ("phi", "dtheta"),
        ("phi", "dphi"),
    ]

    # ---- Create figure ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, (v1, v2) in zip(axes, var_pairs):

        df = pd.DataFrame({
            v1: sol[v1],
            v2: sol[v2]
        })

        sns.kdeplot(data=df, x=v1, y=v2, fill=True, cmap="viridis", levels=40, thresh=0.05, ax=ax, cbar=True)

        ax.set_xlabel(v1)
        ax.set_ylabel(v2)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Phase-Space Density KDE - {time_step}", fontsize=14)
    fig.tight_layout()

    return fig

## 
# ------------------------------ DRIFT ENERRGY GRAPHICS --------------------------
##
def drift_energy(energy:Sequence[float]) -> float:

        return np.abs(energy[-1] - energy[0])

def heatmaps(grid: Dict) -> plt.Figure:

    """
    Generate a heatmap of energy drift for the full (theta, phi) initial condition.

    Paramters:
    ------------
    grid: Sequence
        Nested dictionary grid[theta][phi] = result_dict
    """

    theta_vals = sorted(grid.keys())
    phi_vals = sorted(grid[theta_vals[0]].keys())

    theta, phi = np.meshgrid(phi_vals, theta_vals) # This produces an arrays of shape len(theta_vals) x len(phi_values)
    drift_matrix = np.zeros((len(theta_vals), len(phi_vals)))

   
    for i, th in enumerate(theta_vals):
        for j, ph in enumerate(phi_vals):

            _, _, energy = grid[th][ph]

            energy_drift = drift_energy(energy=energy)

            drift_matrix[i,j] = energy_drift

    fig = plt.figure(figsize = (12,6))

    ax = fig.add_subplot(121, projection="3d")

    # --- 3d plot ------

    surf = ax.plot_surface(theta, phi, drift_matrix, cmap = "viridis", edgecolor = "none")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\phi$")
    ax.set_zlabel(r"$\Delta E$")

    # --- 2d plot ----

    axes = fig.add_subplot(122)
    im = axes.imshow(drift_matrix, extent=[theta_vals[0], theta_vals[-1], phi_vals[0], phi_vals[-1]], origin = "lower", aspect="auto", cmap = "viridis")
    axes.set_xlabel(r"$\phi$")
    axes.set_ylabel(r"$\theta$")
    axes.set_title("Energy drift map")
    fig.colorbar(im, ax = axes)

    plt.tight_layout()

    return fig

def vertical_plane(grid: Dict) -> plt.Figure:
    # --- Extract grid values ---
    theta_vals = sorted(grid.keys())
    phi_vals = sorted(grid[theta_vals[0]].keys())

    # Meshgrid aligned with drift_matrix
    phi, theta = np.meshgrid(phi_vals, theta_vals)
    drift_matrix = np.zeros((len(theta_vals), len(phi_vals)))

    # Compute drift matrix
    for i, th in enumerate(theta_vals):
        for j, ph in enumerate(phi_vals):
            _, _, energy = grid[th][ph]
            drift_matrix[i, j] = drift_energy(energy)

    # --- Pick two random vertical planes ---
    import random
    theta_plane = random.choice(theta_vals)   # vertical plane parallel to phi-axis
    phi_plane   = random.choice(phi_vals)     # vertical plane parallel to theta-axis

    # Extract cross-sections
    # theta = constant -> row in drift_matrix
    drift_theta_cut = drift_matrix[theta_vals.index(theta_plane), :]

    # phi = constant -> column in drift_matrix
    drift_phi_cut = drift_matrix[:, phi_vals.index(phi_plane)]

    # --- Figure layout ---
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(14, 7))
    gs = GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1], figure=fig)

    # --- Main 3D plot ---
    ax1 = fig.add_subplot(gs[:, 0], projection="3d")
    surf = ax1.plot_surface(theta, phi, drift_matrix, cmap="viridis", edgecolor="none", alpha=0.9)

    ax1.set_xlabel(r"$\theta$")
    ax1.set_ylabel(r"$\phi$")
    ax1.set_zlabel(r"$\Delta E$")
    ax1.set_title("Energy Drift Surface with Two Vertical Planes")

    # --- Subplot: θ = θ_plane cut ---
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(phi_vals, drift_theta_cut, color="red")
    ax2.set_xlabel(r"$\phi$")
    ax2.set_ylabel(r"$\Delta E$")
    ax2.set_title(fr"Cut at $\theta = {theta_plane:.3f}$")

    # --- Subplot: φ = φ_plane cut ---
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(theta_vals, drift_phi_cut, color="blue")
    ax3.set_xlabel(r"$\theta$")
    ax3.set_ylabel(r"$\Delta E$")
    ax3.set_title(fr"Cut at $\phi = {phi_plane:.3f}$")

    plt.tight_layout()
    return fig



##
# -----------------------------COMPUTE --------------------------------
##

def compute(data: bool = False, flag:bool = False, stored: bool = False):
    """
    The function is desgined to run the whole error_graphics.py file. The current parameters are:

    Parameters
    -------------
    data: bool
        the paremeter controls whether the store_data.npz file will be use or not. Default = False
    flag: bool
        the parameter controls whether the simulation will be linearized or nor. Default = False
    stored: bool
        the paremeter controls whether the current grpahics are stored or not. Default = False
    
    Friendly reminder
    -----------------
    You can change all parameters you want in manually way. Time, length, mass...
    """

    dt = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
    
    results = run(dt=dt, time= 150, flag = flag)

    print("Starting with figures:")
    runtime = time_compute(dt = dt, results_dict=results)
    conver = convergence(dt = dt, results_dict=results)
    stab = stability(dt = dt, results_dict= results)
    kde = kde_phase_space_subplots(results_dict=results)

    if data:
        print("Starting data files and stored the figures")

        directory = os.getcwd()
        route = os.path.join(directory, "figures")
        os.makedirs(route, exist_ok=True)

        if flag:
            
            grid = upload_data(name = "linearized_stored_data.npz")
            
            heat = heatmaps(grid = grid)
            vertical = vertical_plane(grid=grid)

            conver.savefig(fname= os.path.join(route, "convergence_linearized.png"), dpi = 300, bbox_inches = "tight")
            stab.savefig(fname= os.path.join(route, "stability_linearized.png"), dpi = 300, bbox_inches = "tight")
            kde.savefig(fname= os.path.join(route, "kde_phase_space_linearized.png"), dpi = 300, bbox_inches = "tight")
            heat.savefig(fname= os.path.join(route, "heat_linearized.png"), dpi = 300, bbox_inches = "tight")
            vertical.savefig(fname= os.path.join(route, "vertical_linearized.png"), dpi = 300, bbox_inches = "tight")
            runtime.savefig(fname= os.path.join(route, "time_compute_linearized.png"), dpi = 300, bbox_inches = "tight")

        else:

            grid = upload_data(name = "stored_data.npz")

            heat = heatmaps(grid = grid)
            vertical = vertical_plane(grid = grid)
            
            conver.savefig(fname= os.path.join(route, "convergence.png"), dpi = 300, bbox_inches = "tight")
            stab.savefig(fname= os.path.join(route, "stability.png"), dpi = 300, bbox_inches = "tight")
            kde.savefig(fname= os.path.join(route, "kde_phase_space.png"), dpi = 300, bbox_inches = "tight")
            heat.savefig(fname= os.path.join(route, "heat.png"), dpi = 300, bbox_inches = "tight")
            vertical.savefig(fname= os.path.join(route, "vertical.png"), dpi = 300, bbox_inches = "tight")
            runtime.savefig(fname= os.path.join(route, "time_compute.png"), dpi = 300, bbox_inches = "tight")

    else:   
        plt.plot()


if __name__ == "__main__":

    compute(data = True, stored = True)