"""
MAIN GRPAHICS: double pendulum

Small angles
- Drift energy (E(t) - E(0)) vs intial angle (theta1) vs intial angle (theta2) (heatmap - after a long time)
- Regime summary (position, energies, phase space)
- Verlet vs rk45 vs DOP853 methods

Motion
- Drift energy vs large time
- Regime summary
- Fractal motion
- Lyapunov coefficient
- Poincare sections

Another performance:
- Convergence 
- Stability
- Normal modes animation (small angles): Symmetric mode and Antisymmetric mode
- Trajectory + Poincaré section side-by-side: Shows how the chaotic cloud emerges from the trajectory.

"""

from dataclasses import dataclass
from collections.abc import Sequence
from double_pendulum import DoublePendulum

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Stablish automatically: font sizes, grid visibility, color harmony, spacing
plt.style.use("seaborn-v0_8-paper")
sns.set_theme(context="notebook", style="whitegrid", palette="viridis", font_scale=1.2)


##
#---------------- Parameters and control function -------------------
##
def _as_pair(values: Sequence[float], name: str) -> tuple[float, float]:
    if len(values) != 2:
        raise ValueError(f"{name} must ocntain exactly two values.")
    return float(values[0]), float(values[1])

def _time_grid(duration: float, dt: float) -> np.ndarray:

    if duration <= 0.0:
        raise ValueError("Simulation time must be positive. ")
    if dt <= 0.0:
        raise ValueError("Simulation step must be positive. ")
    
    times = np.arange(0.0, duration, dt, dtype=float)

    #Force exact final time

    times[-1] = duration
    return times

@dataclass
class Params:
    """Physical parameters and initial conditions for the pendulum."""

    g: float = 9.81  # m/s^2

    m1: float = 1.0  # kg
    m2: float = 1.0  # kg
    L1: float = 1.0  # m
    L2: float = 2.0  # m

    q0: tuple[float, float] = (np.deg2rad(10.0), np.deg2rad(10.0))  # rad
    dq0: tuple[float, float] = (0.0, 0.0)  # rad/s

    t: float = 15.0  # s
    dt: float = 0.01  # s

    def __post_init__(self) -> None: #Validate parameters. Default function after dataclass function
        if self.g <= 0.0:
            raise ValueError("Gravity must be positive.")
        if self.m1 <= 0.0 or self.m2 <= 0.0:
            raise ValueError("Masses must be positive.")
        if self.L1 <= 0.0 or self.L2 <= 0.0:
            raise ValueError("Lengths must be positive.")

        self.q0 = _as_pair(self.q0, "q0")
        self.dq0 = _as_pair(self.dq0, "dq0")
        self.times = _time_grid(duration= self.t, dt = self.dt)

###
# --------------------- Main graphics --------------------------
###

def regime_summary(sol:Sequence[float], energy: Sequence[float], position: Sequence[float], colors: Sequence[float], name: str):
    """
    Regime summary -- Position, energies and phase space
    """
    if sol is None:
        raise ValueError("Run the simulation first.")

    theta1, theta2, omega1, omega2 = sol["y"]
    time = sol["t"]
    T, U, Et = energy
    x1, y1, x2, y2 = position

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 7))

    fig.suptitle(f"Double Pendulum Summary — {name}", fontsize=16, weight="bold")
    fig.tight_layout()

    # Position
    ax1.plot(x1, y1, color=colors["mass1"], label="Mass 1")
    ax1.plot(x2, y2, color=colors["mass2"], label="Mass 2")
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlabel(r"$x$ [m]")
    ax1.set_ylabel(r"$y$ [m]")
    ax1.set_title("Position")
    ax1.legend()
    ax1.grid(True, which="both", linestyle="--", alpha=0.4)

    # Energies
    ax2.plot(time, T, color=colors["T"], label="Kinetic")
    ax2.plot(time, U, color=colors["U"], label="Potential")
    ax2.plot(time, Et, color=colors["Et"], label="Total")
    ax2.set_title("Energies")
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Energy [J]")
    ax2.legend()
    ax2.grid(True, which="both", linestyle="--", alpha=0.4)
    ax2.annotate(f"Drift = {Et[-1] - Et[0]:.2e}", xy=(0.05, 0.1), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    # Phase space
    ax3.plot(theta1, omega1, color=colors["mass1"], label="Mass 1")
    ax3.plot(theta2, omega2, color=colors["mass2"], label="Mass 2")
    ax3.set_title("Phase Space")
    ax3.set_xlim(min(theta1.min(), theta2.min()), max(theta1.max(), theta2.max()))
    ax3.set_ylim(min(omega1.min(), omega2.min()), max(omega1.max(), omega2.max()))
    ax3.set_xlabel(r"$\theta$ [rad]")
    ax3.set_ylabel(r"$\omega$ [rad/s]")
    ax3.legend()
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    plt.show()

def initial_angle_error(params, theta1: Sequence[float], theta2: Sequence[float]) -> np.ndarray:

    if len(theta1) or len(theta2) == 0:
        raise ValueError("Theta1 or theta2 must be contains both at least one angle")
    

    angles_theta1 = {angle: {"Drift energy":[]} for angle in theta1}
    angles_theta2 = {angle: {"Drift energy": []} for angle in theta2}
    conbination = {angle: {"Drift energy": []} for angle in theta1}

    for theta in theta1:

        params.q0 = theta, 0
        double =DoublePendulum(params =P, small_angle=True, method = "Verlet")
        sol = double.run()
        energy = double.energies()

        _, _ ,Et = energy
        Energy_drift =np.abs(Et[-1] - Et[0])

        angles_theta1[theta]["Drift energy"].append(Energy_drift)

    for theta in theta2:

        params.q0 = 0, theta
        double =DoublePendulum(params =P, small_angle=True, method = "Verlet")
        sol = double.run()
        energy = double.energies()

        _, _ ,Et = energy
        Energy_drift =np.abs(Et[-1] - Et[0])

        angles_theta2[theta]["Drift energy"].append(Energy_drift)

    

P = Params()
double = DoublePendulum(params= P, small_angle=True, method="Verlet")
sol = double.run()
postion = double.transform()
energy = double.energies()

cmap = plt.colormaps["viridis"]
color = cmap(np.linspace(0, 1, 3))
colors = {
    "mass1": cmap(0.2),
    "mass2": cmap(0.8),
    "T": "#1f77b4",
    "U": "#ff7f0e",
    "Et": "#2ca02c",
}

regime_summary(sol = sol, energy=energy, position=postion, colors = colors, name = "Verlet")


"""
class SmallAngles:
    def __init__(self, Params: Sequence[float], Instance: float):

        self.Params = Params
        self.Instance = Instance
    def run():

    def phase_space():
    def drift_energy_surface():
    def method_error():
"""
