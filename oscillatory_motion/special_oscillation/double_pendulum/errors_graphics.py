"""
MAIN ERRORS GRAPHICS: double pendulum

Linearized equation:
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Drift energy vs large time ----- Done

Normal equation:
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Drift energy vs large time ----- Done

Numerical methods:
- Convergence 
- Stability
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

    theta1_0: float = np.deg2rad(145)
    theta2_0: float = np.deg2rad(45)

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

params = Params()
sim = DoublePendulumSimulator(params = params)


def drift_energy_comparison(results: Sequence[float]):
    
    fig = plt.figure(figsize=(10, 6))
    #Energy comparisons
    colors = ['blue', 'green', 'red', 'orange']
    
    for i, (name, (result, _)) in enumerate(results.items()):
        plt.plot(result.t, result.total - result.total[0], 
                 label=name, color=colors[i], linewidth=1.5)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Energy Drift (J)')
    plt.title('Total Energy Drift Over Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('symlog')  # Symmetric log to show small drifts

    return fig


results_dict = {}

methods_to_test = [
        ("RK45 (adaptive explicit)", "RK45", False),
        ("DOP853 (8th order)", IntegrationMethod.DOP853, False),
        ("Radau (implicit)", IntegrationMethod.RADAU, False),
        ("BDF (backward diff)", IntegrationMethod.BDF, False),
    ]

#Run the all methods.
for name, method, _ in methods_to_test:
        print(f"\nRunning with {name}...")
        result = sim.run(method=method, linearized=False)
        analysis = sim.energy_analysis(result)
        results_dict[name] = (result, analysis)
        
        print(f"  Initial Energy: {analysis['initial_energy']:.6f} J")
        print(f"  Final Energy:  {analysis['final_energy']:.6f} J")
        print(f"  Energy Drift: {analysis['energy_drift']:.6e} J")
        print(f"  Relative Error: {analysis['relative_error_ppm']:.2f} ppm")
    
fig3 = drift_energy_comparison(results = results_dict)