"""
MAIN ERRORS GRAPHICS: double pendulum

Linearized equation:
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Drift energy vs large time
- Time

Normal equation:
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Drift energy vs large time 
- Time 

Numerical methods:
- Convergence 
- Stability
- Phase space densitiy
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
    theta2_0: float = np.deg2rad(0)
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
    Drift_matrix = np.zeros_like(Theta1, dtype=float)
    
    # Populate matrix
    for i, t1 in enumerate(theta1_vals):
        for j, t2 in enumerate(theta2_vals):
            Drift_matrix[j, i] = grid[t1][t2]['drift']

    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Use log scale if drifts vary by orders of magnitude
    # Clamp small values to avoid log(0)
    data_to_plot = np.maximum(Drift_matrix, 1e-15)
    
    sns.heatmap(data_to_plot, 
                xticklabels=np.round(theta1_vals, 2), 
                yticklabels=np.round(theta2_vals, 2),
                cmap="viridis", ax=ax, 
                cbar_kws={'label': 'Energy Drift [J]'},
                norm=plt.matplotlib.colors.LogNorm())
    
    ax.set_title("Energy Drift Heatmap (log scale)")
    ax.set_xlabel(r"Initial $\theta_1$ [rad]")
    ax.set_ylabel(r"Initial $\theta_2$ [rad]")
    
    # Rotate labels for better readability
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=0)
    
    plt.tight_layout()
    return fig


def drift_energy_comparison(dt: Sequence[float], results_dict: Dict) -> plt.Figure:
    """
    Plot the drift in total energy over time for different methods.
    """
    fig = plt.figure(figsize=(10, 6))
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#ff7f0e'] # Standard matplotlib colors
    
    for i, (method_name, dt_dict) in enumerate(results_dict.items()):
        color = colors[i % len(colors)]

        for time_step in dt:
            if time_step not in dt_dict:
                continue

            result, analysis = dt_dict[time_step]
            
            # Access trajectory: assume result is an object with .t and .total attributes
            # Adjust based on your specific DoublePendulumSimulator output
            t = result.t
            E = result.total
            drift = np.abs(E - E[0])

            plt.plot(t, drift, label=f"{method_name} (dt={time_step})", color=color, linewidth=1.5, alpha=0.8)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Energy Drift |E(t) - E(0)| [J]')
    plt.title('Total Energy Drift Over Time')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.yscale('symlog') # Symmetric log to show small drifts

    return fig


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

        plt.plot(dts_sorted, times_sorted, 'o-', label=method_name, color=color, linewidth=2, markersize=6)
    
    plt.xlabel('Time Step ($dt$)')
    plt.ylabel('Run Time (s)')
    plt.title('Computational Cost vs Time Step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    plt.yscale('log')
    
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
                "error": error,
                "method": method_name
            })

    df = pd.DataFrame(rows)

    fig = plt.figure(figsize=(10, 6))
    sns.lineplot(data=df, x="dt", y="error", hue="method", marker="o", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Time Step $dt$")
    plt.ylabel("Energy Drift (Error)")
    plt.title("Convergence Study: Error vs Time Step")
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


def kde_phase_space_subplots(results_dict: Dict, var1: str="theta1", var2: str="omega1") -> plt.Figure:
    """
    Phase-space density plots (KDE) for each integration method.
    """
    n_methods = len(results_dict)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, (method_name, dt_dict) in zip(axes, results_dict.items()):
        # Choose the smallest dt (best resolution) for the cleanest phase plot
        time_step = min(dt_dict.keys())
        
        result, analysis = dt_dict[time_step]
        
        sol = result.y
        
        if hasattr(sol, var1):
            x = getattr(sol, var1)
            y = getattr(sol, var2)
        elif hasattr(sol, '__getitem__'):
            # Assuming indexed access: y[0] is theta1, y[1] is theta2 etc.
            # Map indices to variable names
            idx_map = {"theta1": 0, "theta2": 1, "omega1": 2, "omega2": 3}
            x = sol[idx_map[var1]]
            y = sol[idx_map[var2]]
        else:
            print(f"Cannot extract data for {var1}, {var2} from result object.")
            continue

        # Create DataFrame for Seaborn KDE
        df = pd.DataFrame({
            var1: x,
            var2: y
        })

        # Plot KDE
        sns.kdeplot(data=df, x=var1, y=var2, fill=True, cmap="viridis", levels=40, thresh=0.05, ax=ax)
        
        ax.set_title(f"Phase Space: {method_name} (dt={time_step})")
        ax.set_xlabel(var1)
        ax.set_ylabel(var2)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Phase-Space Density: {var1} vs {var2}", fontsize=16)
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
        ("RK45 (adaptive explicit)", IntegrationMethod.RK45, False),
        ("DOP853 (8th order)", IntegrationMethod.DOP853, False),
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
    
    # 1. Run Convergence & Energy Drift Analysis
    dt_values = [0.5, 0.1, 0.05, 0.01]
    
    print("Running simulations...")
    results = compute(dt=dt_values, time=50, flag=False)
    
    # Plotting Numerical Methods Results
    print("Plotting Convergence")
    fig_conv = convergence(dt=dt_values, results_dict=results)
    plt.show()
    print("Plotting Drift energy comparison")
    fig_drift = drift_energy_comparison(dt=dt_values, results_dict=results)
    plt.show()
    print("Plotting time compute")
    fig_time = time_compute(dt=dt_values, results_dict=results)
    plt.show()
    print("Plotting phase spacev theta1 vs omega1")
    fig_phase = kde_phase_space_subplots(results_dict=results, var1="theta1", var2="omega1")
    plt.show()
    print("Plotting phase spacev theta2 vs omega2")
    fig_phase = kde_phase_space_subplots(results_dict=results, var1="theta2", var2="omega2")
    plt.show()
    
    # 2. Load Initial Condition Data
    try:
        theta1_s, theta2_s, grid_s = upload_data("compute_initial_angles.npz")
        
        # Initial Angle Drift Plots
        fig_init = initial_angle_drift_energy(theta1_s, theta2_s)
        plt.show()
        
        # Heatmaps
        fig_heat = heatmaps(grid_s)
        plt.show()
        
    except FileNotFoundError as e:
        print(e)
        print("Skipping heatmaps. Run compute_data.py to generate data.")