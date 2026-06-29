"""
MAIN ERRORS GRAPHICS: double_pendulum.py 

This code generate the following graphics such as linearizad equation as normal equation

- Drift energy (E(t) - E(0)) vs intial angle (theta1)
- Drift energy (E(t) - E(0)) vs intial angle (theta2)
- Heatmap 
- Drift energy vs large time
- Time

Numerical methods:
- Convergence 
- Stability
- Phase space densitiy

I fervently recommend you run the code with different parameters
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
from double_pendulum import DoublePendulumSimulator
from dataclasses import dataclass
from typing import Sequence, Dict, Any, Tuple, List
from enum import Enum, auto

##
# ------------------- PREPARATIONS ------------------
##
class IntegrationMethod(Enum):
    """Supported integration methods with their properties."""
    RK45 = auto()          # 4th order explicit - not symplectic, energy drifts
    DOP853 = auto()        # 8th order explicit - better than RK45
    RADAU = auto()         # Implicit Radau - BEST for energy conservation
    BDF = auto()           # Backward Differentiation - good for stiff systems

@dataclass
class Params:
    """Physical parameters and initial conditions for the pendulum."""
    g: float = 9.81  # m/s^2
    m1: float = 1.0  # kg
    m2: float = 1.0  # kg
    L1: float = 1.0  # m
    L2: float = 2.0  # m

    theta1_0: float = np.deg2rad(45)
    theta2_0: float = np.deg2rad(45)
    omega1_0: float = 0.0
    omega2_0: float = 0.0

    t_max: float = 15.0  # s
    dt: float = 1e-3  # s
    rtol: float = 1e-10
    atol: float = 1e-12

    def __post_init__(self) -> None:
        if self.g <= 0.0:
            raise ValueError("Gravity must be positive.")
        if self.m1 <= 0.0 or self.m2 <= 0.0:
            raise ValueError("Masses must be positive.")
        if self.L1 <= 0.0 or self.L2 <= 0.0:
            raise ValueError("Lengths must be positive.")

##
# ------------------------- UPLOAD DATA --------------------------
##
def upload_data(filename: str = "compute_intial_angles.npz") -> Tuple[Dict, Dict, Dict]:
    """
    Load precomputed simulation data from compute_data.py.
    
    Returns
    -------
    Tuple[Dict, Dict, Dict]
        (theta1_scan, theta2_scan, grid)
        Each dictionarycontains keys mapped to result dictionaries with 'drift', 'theta1', etc.
    """
    route = os.path.join(os.getcwd(), filename)
    
    if not os.path.isfile(route):
        raise FileNotFoundError(f"Data file not found: {route}. Run compute_data.py first.")

    # Load compressed numpy dictionary
    data = np.load(route, allow_pickle=True)
    
    # Extract the specific dictionary items
    # These correspond to the kwargs saved in compute_data.py
    theta1_scan = data["theta1_scan"].item()
    theta2_scan = data["theta2_scan"].item()
    grid = data["grid"].item()
    
    return theta1_scan, theta2_scan, grid

##
# ----------------------- STUDY OF ENERGY DRIFT ----------------------
##

def initial_angle_drift_energy(theta1_scan: Dict, theta2_scan: Dict) -> plt.Figure:
    """
    Analyze how energy drift depends on initial angles.
    
    Parameters
    ----------
    theta1_scan : Dict
        Results of varying theta1 (theta2=0). Keys are angles (floats), values are result dicts.
    theta2_scan : Dict
        Results of varying theta2 (theta1=0). Keys are angles (floats), values are result dicts.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Process Theta 1 Scan
    # Extract sorted angles and corresponding drifts
    angles1 = sorted(theta1_scan.keys())
    drifts1 = [theta1_scan[ang]['drift'] for ang in angles1]
    
    axes[0].plot(angles1, drifts1, 'o-', color='blue', linewidth=2, markersize=6)
    axes[0].set_title(r"Energy Drift vs Initial $\theta_1$ ($\theta_2=0$)")
    axes[0].set_xlabel(r"Initial Angle $\theta_1$ [rad]")
    axes[0].set_ylabel("Maximum Energy Drift [J]")
    axes[0].grid(True, alpha=0.3)
    
    # Process Theta 2 Scan
    angles2 = sorted(theta2_scan.keys())
    drifts2 = [theta2_scan[ang]['drift'] for ang in angles2]
    
    axes[1].plot(angles2, drifts2, 'o-', color='green', linewidth=2, markersize=6)
    axes[1].set_title(r"Energy Drift vs Initial $\theta_2$ ($\theta_1=0$)")
    axes[1].set_xlabel(r"Initial Angle $\theta_2$ [rad]")
    axes[1].set_ylabel("Maximum Energy Drift [J]")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def heatmaps(grid: Dict) -> plt.Figure:
    """
    Generate heatmaps of energy drift for the full (theta1, theta2) parameter space.
    
    Parameters
    ----------
    grid : Dict
        Nested dictionary grid[theta1][theta2] = result_dict
    """
    # Extract unique sorted angles
    theta1_vals = sorted(grid.keys())
    theta2_vals = sorted(grid[theta1_vals[0]].keys())
    
    # Create meshgrid matrices
    Theta1, Theta2 = np.meshgrid(theta1_vals, theta2_vals)
    Drift_matrix = np.zeros((len(theta1_vals), len(theta2_vals)))
    
    # Populate matrix
    for i, t1 in enumerate(theta1_vals):
        for j, t2 in enumerate(theta2_vals):
            Drift_matrix[i, j] = grid[t1][t2]['drift']

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121, projection = "3d") 
    # ------- 3D plot --------
    surf = ax.plot_surface(Theta1, Theta2, Drift_matrix, cmap = "viridis", edgecolor = "none")
    ax.set_xlabel(r'$\theta_{0_1}$')
    ax.set_ylabel(r'$\theta_{0_2}$')
    ax.set_title(f'Energy Drift Surface')

    # ---- 2D HEATMAP ----
    ax2 = fig.add_subplot(122)
    im = ax2.imshow(Drift_matrix, extent=[theta2_vals[0], theta2_vals[-1], theta1_vals[0], theta1_vals[-1]] ,origin='lower', aspect='auto', cmap='viridis')
    ax2.set_xlabel(r'$\theta_{0_1}$')
    ax2.set_ylabel(r'$\theta_{0_2}$')
    ax2.set_title(f'Energy Drift Map')
    fig.colorbar(im, ax=ax2)

    plt.tight_layout()
    return fig


def drift_energy_comparison(dt: Sequence[float], results_dict: Dict, method_name: str = "DOP853 (8th order)") -> Tuple[plt.Figure, plt.Figure]:
    """
    Plot the drift in total energy over time.

    fig1:
        Same method, different time steps.
    fig2:
        Same time step (min dt) across different methods.
    """

    colors = ['#1f77b4', '#2ca02c', '#d62728', "#ff0edf"]

    # ============================================================
    # FIGURE 1 — SAME METHOD, DIFFERENT dt
    # ============================================================
    fig1, ax1 = plt.subplots(figsize=(10, 6))

    for i, dt_value in enumerate(dt):

        if dt_value not in results_dict[method_name]:
            continue

        result, analysis = results_dict[method_name][dt_value]

        t = result.t
        E = result.total
        drift = np.abs(E - E[0])

        ax1.plot(t, drift, label=f"dt={dt_value}", color=colors[i % len(colors)], linewidth=1.5,alpha=0.8)

    ax1.set_title(f"Total Energy Drift Over Time — {method_name}")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel(r"Energy Drift $|E(t) - E(0)|$ [J]")
    ax1.set_xlim([0,15])
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale("symlog")

    # ============================================================
    # FIGURE 2 — SAME dt, DIFFERENT METHODS
    # ============================================================
    fig2, ax2 = plt.subplots(figsize=(10, 6))

    # choose the smallest dt available for the reference method
    min_step = min(results_dict[method_name].keys())

    for i, (method, dt_dict) in enumerate(results_dict.items()):

        if min_step not in dt_dict:
            continue

        result, analysis = dt_dict[min_step]

        t = result.t
        E = result.total
        drift = np.abs(E - E[0])

        ax2.plot(t, drift, label=f"{method}", color=colors[i % len(colors)], linewidth=1.5, alpha=0.8)

    ax2.set_title(f"Total Energy Drift Across Methods — dt={min_step}")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel(r"Energy Drift $|E(t) - E(0)|$ [J]")
    ax2.set_xlim([0,15])
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale("symlog")

    return fig1, fig2



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
            
            result, analysis = dt_dict[time_step]
            
            # Access computation time if available, otherwise estimate via array length
            # Adjust 'result.total_time' or similar based on your simulator's output
            if hasattr(result, 'total_time'):
                run_time = result.total_time
            else:
                # Fallback: simply mark based on array size * dt
                run_time = len(result.t) * time_step

            times.append(run_time)
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

##
# ---------------------------- NUMERICAL STUDY METHODS ------------------------
##

def convergence(dt: Sequence[float], results_dict: Dict) -> plt.Figure:
    """
    Study numerical convergence: Error vs Time Step.
    """
    rows = []

    for method_name, dt_dict in results_dict.items():
        for time_step in dt:
            if time_step not in dt_dict:
                continue

            result, analysis = dt_dict[time_step]
            error = analysis["energy_drift"] 

            rows.append({
                "dt": time_step,
                "error": np.abs(error),
                "method": method_name
            })

    df = pd.DataFrame(rows)

    fig = plt.figure(figsize=(10, 6))
    sns.lineplot(data=df, x="dt", y="error", hue="method", marker="o", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Time Step $dt$")
    plt.ylabel("Energy Drift (Error)")
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
                
            result, analysis = dt_dict[time_step]
            
            # Check trajectory for NaNs or Infs
            E = result.total
            has_nan = np.any(np.isnan(E))
            has_inf = np.any(np.isinf(E))
            
            # Determine stability label
            status = "Unstable" if (has_nan or has_inf) else "Stable"
            
            rows.append({
                "dt": time_step,
                "method": method_name,
                "status": status,
                "final_energy": E[-1] if not has_nan and not has_inf else np.nan
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
    time_step = min(results_dict["DOP853 (8th order)"].keys())
    result, analysis = results_dict["DOP853 (8th order)"][time_step]

    # ---- Extract solution arrays ----
    # Case 1: Your custom result object with attributes
    if hasattr(result, "theta1"):
        sol = {
            "theta1": result.theta1,
            "theta2": result.theta2,
            "omega1": result.omega1,
            "omega2": result.omega2,
        }

    # Case 2: SciPy OdeResult (result.y is array)
    elif hasattr(result, "y"):
        idx_map = {"theta1": 0, "theta2": 1, "omega1": 2, "omega2": 3}
        sol = {name: result.y[idx] for name, idx in idx_map.items()}

    else:
        raise ValueError("Result object has no recognizable structure.")

    # ---- Variables to plot ----
    var_pairs = [
        ("theta1", "omega1"),
        ("theta1", "omega2"),
        ("theta2", "omega1"),
        ("theta2", "omega2"),
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
# -------------------------- Compute each instance -----------------------
##
def compute(dt: Sequence[float], time: int = 150, flag: bool = False,) -> Dict:
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

    methods_to_test = [
        ("DOP853 (8th order)", IntegrationMethod.DOP853, False),
        ("RK45 (adaptive explicit)", IntegrationMethod.RK45, False),
        ("Radau (implicit)", IntegrationMethod.RADAU, False),
        ("BDF (backward diff)", IntegrationMethod.BDF, False),
    ]

    #Large time
    params.t_max = time

    #Instance
    sim = DoublePendulumSimulator(params = params)

    results_dict = {}
    #Run the all methods.
    for name, method, _ in methods_to_test:

        results_dict[name] = {}
        for time_step in dt:

            params.dt = time_step

            print(f"\nRunning with {name} with {time_step}...")
            result = sim.run(method=method, linearized=flag)
            analysis = sim.energy_analysis(result)

            results_dict[name][time_step] = (result, analysis)
        
            print(f"  Initial Energy: {analysis['initial_energy']:.6f} J")
            print(f"  Final Energy:  {analysis['final_energy']:.6f} J")
            print(f"  Energy Drift: {analysis['energy_drift']:.6e} J")
            print(f"  Relative Error: {analysis['relative_error_ppm']:.2f} ppm")
    
    return results_dict

# --- Main Execution Block ---

if __name__ == "__main__":
    
    stored = True
    # 1. Run Convergence & Energy Drift Analysis
    dt_values = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
    
    print("Running simulations...")
    results = compute(dt=dt_values, time=50, flag=True)# Change the flag to compute the non-linearized system
    
    # Plotting Numerical Methods Results
    print("Plotting Convergence")
    fig_conv = convergence(dt=dt_values, results_dict=results)
    
    print("Plotting stability")
    fig_stab = stability(dt= dt_values, results_dict= results) # mirar
        
    print("Plotting Drift energy comparison")
    fig_drift_dt, fig_drift_method = drift_energy_comparison(dt=[1, 0.1, 0.01, 0.001], results_dict=results)
    
    print("Plotting time compute")
    fig_time = time_compute(dt=dt_values, results_dict=results)

    print("Plotting phase space")
    fig_phase1 = kde_phase_space_subplots(results_dict=results)
    
    # 2. Load Initial Condition Data
    """
    try:
        theta1_s, theta2_s, grid_s = upload_data("compute_intial_angles.npz")
        
        # Initial Angle Drift Plots
        fig_init = initial_angle_drift_energy(theta1_s, theta2_s)
        
        # Heatmaps
        fig_heat = heatmaps(grid_s)
        
    except FileNotFoundError as e:
        print(e)
        print("Skipping heatmaps. Run compute_data.py to generate data.")
    """
    if stored:
        directory = os.getcwd()
        route = os.path.join(directory, "figures")
        os.makedirs(route, exist_ok=True)

        fig_conv.savefig(os.path.join(route, "convergence_linearized.png"), dpi = 300, bbox_inches = "tight")
        fig_drift_method.savefig(os.path.join(route, "drift_energy_method_linearized.png"), dpi = 300, bbox_inches = "tight")
        fig_drift_dt.savefig(os.path.join(route, "drift_energy_dt_linearized.png"), dpi = 300, bbox_inches = "tight")
        fig_phase1.savefig(os.path.join(route, "phase_density_linearized.png"), dpi = 300, bbox_inches = "tight")
        fig_time.savefig(os.path.join(route, "runtime_linearized.png"), dpi = 300, bbox_inches = "tight")
        fig_stab.savefig(os.path.join(route, "stability_linearized.png"), dpi = 300, bbox_inches = "tight")
        #fig_init.savefig(os.path.join(route, "initial_condition.png"), dpi = 300, bbox_inches = "tight")
        #fig_heat.savefig(os.path.join(route, "heatmap.png"), dpi = 300, bbox_inches = "tight")