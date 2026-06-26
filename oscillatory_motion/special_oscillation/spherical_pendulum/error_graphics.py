"""
Spherical pendulum - main graphics:

Numerical graphics can be run by both: linearized and normal simulation.

- Comparison method.
- Energy drift colormap.
- Runtime. ---> Done
- Convergence/stability. --> Done
- Density kde ---> done

Special grahics:
- Invariant-torus reconstruction. Use delay embedding or Fourier decomposition to visualize the torus in 3D.
- Frequency-map analysis (Laskar). Plot frequency drift to detect weak chaos. (MIRA EN FAVORITOS DE GOOGLE)
- Action-angle coordinate plots. If you compute approximate actions, J_theta,J_phi, plot trajectories in action space.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import pandas as pd

from spherical_pendulum import Spherical_Pendulum
from dataclasses import dataclass
from typing import Tuple, Sequence, Any, Dict
from enum import Enum, auto

def _as_pair(values: Sequence[float], name: str) -> tuple[float, float]:
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
# ------------------------- Compute -------------------------------
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

def upload_data(name:str = "loaded_data.npz", directory: str = os.getcwd()) -> Tuple[Any, Any, Any]:

    path = os.path.join(directory, name)

    if os.path.exists(path):
        
        data = np.load(path, allow_pickle= True)

        theta = data["theta_scan"].item()
        phi = data["phi_scan"].item()
        grid =data["grid"].item()

        return theta, phi, grid

    else:
        raise FileNotFoundError("Compute first the load_data.py file")
    

##
# ------------------ Main Errors Graphics ---------
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
    # Case 1: Your custom result object with attributes
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




results = compute(dt = [1, 0.1, 0.01])

fig = time_compute(dt = [1, 0.1, 0.01], results_dict=results)
fig1 = convergence(dt = [1, 0.1, 0.01], results_dict=results)
fig2 = stability(dt = [1, 0.1, 0.01], results_dict= results)
#fig3 = kde_phase_space_subplots(results_dict=results)
plt.show()