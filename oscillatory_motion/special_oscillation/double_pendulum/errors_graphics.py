"""
MAIN ERRORS GRAPHICS: double pendulum

Linearized equation:
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Drift energy vs large time ----- Done
- Time ----- Done

Normal equation:
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Drift energy vs large time ----- Done
- Time ----- Done

Numerical methods:
- Convergence 
- Stability
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

from dataclasses import dataclass
from typing import Sequence, Dict, Any
from double_pendulum import DoublePendulumSimulator
from enum import Enum, auto   


class IntegrationMethod(Enum):
    """Supported integration methods with their properties."""
    RK45 = auto()           # 4th order explicit - not symplectic, energy drifts
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

    rtol = 1e-10
    atol = 1e-12

    def __post_init__(self) -> None: #Validate parameters. Default function after dataclass function
        if self.g <= 0.0:
            raise ValueError("Gravity must be positive.")
        if self.m1 <= 0.0 or self.m2 <= 0.0:
            raise ValueError("Masses must be positive.")
        if self.L1 <= 0.0 or self.L2 <= 0.0:
            raise ValueError("Lengths must be positive.")


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

def drift_energy_comparison(dt: Sequence[float], results_dict: Sequence[float]):

    fig = plt.figure(figsize=(10, 6))

    colors = ['blue', 'green', 'red', 'orange']
    
    for i, (method_name, dt_dict) in enumerate(results_dict.items()):
        color = colors[i % len(colors)]

        for time_step in dt:
            if time_step not in dt_dict:
                continue  # skip missing dt values

            result, analysis = dt_dict[time_step]

            # Energy drift over time
            E = result.total
            drift = np.abs(E - E[0])

            plt.plot(result.t, drift, label=f"{method_name} (dt={time_step})", color=color, linewidth=1.5, alpha=0.8)
    
    plt.xlabel('Time')
    plt.ylabel('Energy Drift')
    plt.title('Total Energy Drift Over Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('symlog')  # Symmetric log to show small drifts

    return fig

def time_compute(dt: Sequence[float], results_dict: Sequence[float]):

    fig = plt.plot(figsize = (10, 6))

    colors = ['blue', 'green', 'red', 'orange']

    for i, (method_name, dt_dict) in enumerate(results_dict.items()):

        color = colors[i % len(colors)]

        for time_step in dt:
            if time_step not in dt_dict:
                continue  # skip missing dt values

            result, analysis = dt_dict[time_step]

            # Energy drift over time
            Time = result.total_time

            plt.plot(time_step, Time, label=f"{method_name} (dt={time_step})", color=color, linewidth=1.5, alpha=0.8)
    
    plt.xlabel('Time step')
    plt.ylabel('Run time')
    plt.title('Run time vs time step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('symlog')
    return fig

##
# -----------------Heatmaps and instability through initial condition -------------------
##

def upload_data(filename: str = "compute_intial_angles.npz "):
    
    route = os.path.join(os.getcwd(), filename)

    if os.path.isfile(route):

        data = np.load(route, allow_pickle= True)

        theta_angle_1 = data["theta1_scan"]
        theta_angle_2 = data["theta2_scan"]
    else:

        print("Compute_data.py must be run first")

    return theta_angle_1, theta_angle_2

def initial_angle_drift_energy(theta1_scan: Dict, theta2_scan: Dict) -> plt.Figure:
    """
    Description
    -----------
    Compute how the drift energy change over different initial theta condition

    Parameters
    ----------
    theta1_scan : Dict,
        Dictionary with these entries: "theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []

    theta2_scan : Dict,
        Dictionary with these entries: "theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []

    Request
    ---------
    Compute first the compute_data.py file. Remember: The greater the number of thetas, the greater the accuracy
    """

    fig = plt.figure(figsize=(10, 6))

    plt.subplot(121)
    plt.plot(theta1_scan, theta1_scan.Drift, label = "Theta 1")
    plt.title("Initial theta 1 vs drift energy")
    plt.xlabel(r"Theta 1 [rad]")
    plt.ylabel(r"Energy drift")

    plt.subplot(122)
    plt.plot(theta2_scan, theta2_scan.Drift, label = "Theta 2")
    plt.xlabel(r"Theta 2 [rad]")
    plt.ylabel(r"Energy drift")
    plt.title("Initial theta 2 vs drift energy")

    plt.tight_layout()

    return fig

def heatmaps(theta_scan: Dict) -> plt.Figure:
    pass

##
# ----------------------- Numerical analysis -----------------
##

results_dict = compute(dt = [0.1, 0.01, 0.005, 0.001])
fig1 = drift_energy_comparison(dt = [0.1, 0.01, 0.001], results_dict= results_dict)
fig2 = time_compute(dt = [0.1, 0.01, 0.005, 0.001], results_dict= results_dict)
plt.show(fig1, fig2)