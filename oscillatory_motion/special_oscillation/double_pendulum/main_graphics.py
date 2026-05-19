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
from typing import Sequence, Dict, Any
from double_pendulum import DoublePendulum

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation
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

    q0: tuple[float, float] = (np.deg2rad(60.0), np.deg2rad(0.0))  # rad
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

def regime_summary(sol:Sequence[float], energy: Sequence[float], position: Sequence[float], colors: Sequence[float], name: str)-> plt.figure:
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
    plt.tight_layout()
    return fig


def double_pendulum_animation(sol:Sequence[float], energy:Sequence[float], position:Sequence[float], colors:Sequence[float], name: str)-> plt.figure:

    theta1, theta2, omega1, omega2 = sol["y"]
    times = sol["t"]
    x1, y1, x2, y2 = position
    T, U, Et = energy

    if len(times) != len(x1) or len(x1) != len(T):
        raise TypeError("Something goes wrong. The length between energy, position and sol don't fix it.")
    
    fig, axes = plt.subplots(1, 3, figsize = (14, 10), tight_layout = True)
    
    # --- MOTION PANEL ---
    ax_m = axes[0]
    rod1, = ax_m.plot([], [], "k-", lw=2)
    rod2, = ax_m.plot([], [], "k-", lw=2)
    m1, = ax_m.plot([], [], "o", color=colors["mass1"], markersize=10)
    m2, = ax_m.plot([], [], "o", color=colors["mass2"], markersize=10)
    ax_m.set_title("Double pendulum motion")
    ax_m.set_xlim(min(x2)-1, max(x2)+1)
    ax_m.set_ylim(min(y2)-1, 1)
    ax_m.set_xlabel("x [m]")
    ax_m.set_ylabel("y [m]")
    ax_m.axhline(0, color="black", lw=0.5)

    # --- ENERGY PANEL ---
    ax_e = axes[1]
    lineT, = ax_e.plot([], [], color=colors["T"], lw=2, label="T")
    lineU, = ax_e.plot([], [], color=colors["U"], lw=2, label="U")
    lineEt, = ax_e.plot([], [], color=colors["Et"], lw=2, label="Et")
    ax_e.set_title("Energies")
    ax_e.set_xlim(times[0], times[-1])
    ax_e.set_ylim(min(U)-1, max(T)+1)
    ax_e.set_xlabel(r"t [s]")
    ax_e.set_ylabel(r"E [J]")
    ax_e.legend()

    # --- PHASE SPACE PANEL ---
    ax_p = axes[2]
    ps1, = ax_p.plot([], [], "--", color=colors["mass1"], lw=2)
    ps2, = ax_p.plot([], [], "--", color=colors["mass2"], lw=2)
    dot1, = ax_p.plot([], [], "o", color=colors["mass1"])
    dot2, = ax_p.plot([], [], "o", color=colors["mass2"])

    combine_theta = theta1 + theta2
    combine_omega = omega1 + omega2
    max_x = max(combine_theta)
    max_y = max(combine_omega)
    min_x = min(combine_theta)
    min_y = min(combine_omega)
    ax_p.set_title(r"Phase space $\theta$ vs $\omega$")
    ax_p.set_xlim(min_x, max_x)
    ax_p.set_ylim(min_y, max_y)
    ax_p.set_xlabel(r"$\theta$ [rad]")
    ax_p.set_ylabel(r"$\omega$ [rad/s]")

    # --- UPDATE FUNCTION ---
    def update(i):
        # Motion
        rod1.set_data([0, x1[i]], [0, y1[i]])
        rod2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        m1.set_data([x1[i]], [y1[i]])
        m2.set_data([x2[i]], [y2[i]])

        # Energies
        lineT.set_data([times[:i]], [T[:i]])
        lineU.set_data([times[:i]], [U[:i]])
        lineEt.set_data([times[:i]], [Et[:i]])

        # Phase space
        ps1.set_data([theta1[:i]], [omega1[:i]])
        ps2.set_data([theta2[:i]], [omega2[:i]])
        dot1.set_data([theta1[i]], [omega1[i]])
        dot2.set_data([theta2[i]], [omega2[i]])

        return rod1, rod2, m1, m2, lineT, lineU, lineEt, ps1, ps2, dot1, dot2

    frame_step = 5
    anim = FuncAnimation(fig, update, frames = np.arange(0, len(times), frame_step), interval = 50, blit = False, repeat = False)
    return anim

def compute_initial_angles(params, theta1: Sequence[float], theta2: Sequence[float], filename: str = (
        r"C:\Users\JORGE\Desktop\Programas\Python\physics_simulation" 
        r"\oscillatory_motion\special_oscillation\double_pendulum\initial_small_angle.npz")
    ) -> Dict[str, Any]:
    """
    Compute energy drift for:
    - varying theta1 with theta2 = 0
    - varying theta2 with theta1 = 0
    - full grid (theta1, theta2)

    Returns a dictionary with three entries:
        "theta1_scan", "theta2_scan", "grid"
    """
    #Validation
    if len(theta1) == 0 or len(theta2) == 0:
        raise ValueError("theta1 and theta2 must each contain at least one angle")

    #Define dictionaries
    theta1_scan = {th: None for th in theta1}
    theta2_scan = {th: None for th in theta2}
    grid = {th1: {th2: None for th2 in theta2} for th1 in theta1}

    
    def _compute_energy_drift(q0):
        params.q0 = q0
        double = DoublePendulum(params=params, small_angle=True, method="Verlet")
        double.run()
        _, _, Et = double.energies()
        return float(abs(Et[-1] - Et[0]))

    # Scan theta1 (theta2 = 0)
    for th in theta1:
        theta1_scan[th] = _compute_energy_drift((th, 0.0))

    #Scan theta2 (theta1 = 0)
    for th in theta2:
        theta2_scan[th] = _compute_energy_drift((0.0, th))

    # Full grid (theta1, theta2)
    total = len(theta1)
    bar_len = 20

    for i, th1 in enumerate(theta1):
        progress = (i + 1) / total
        filled = int(progress * bar_len)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(f"[{bar}]  {progress*100:5.1f}%   θ₁ = {th1:.4f}", end="\r", flush=True)

        for th2 in theta2:
            grid[th1][th2] = _compute_energy_drift((th1, th2))

    print()  # newline after progress bar

    # -----------------------------
    # 6. Save results
    # -----------------------------
    np.savez(filename, theta1_scan=theta1_scan, theta2_scan=theta2_scan, grid=grid)

    return {"theta1_scan": theta1_scan, "theta2_scan": theta2_scan, "grid": grid}

    
P = Params()
double = DoublePendulum(params= P, small_angle=True, method="Rk45")
sol = double.run()
position = double.transform()
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

#regime_summary(sol = sol, energy=energy, position=postion, colors = colors, name = "Verlet")
fig = double_pendulum_animation(sol = sol, energy=energy, position=position, colors= colors, name = "Rk45")
plt.show()

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
