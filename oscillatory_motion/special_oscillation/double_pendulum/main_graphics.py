"""
MAIN GRPAHICS: double pendulum

Linearized equation

- Regime summary (position, energies, phase space) ---Done
- Animation ---- done
- Fractal 

Normal equation

- Regime summary -------- Done
- Fractal 
- Lyapunov coefficient ----- Done
- Poincare sections

Another performance:
- Times
- Normal modes animation (small angles): Symmetric mode and Antisymmetric mode
- Trajectory + Poincaré section side-by-side: Shows how the chaotic cloud emerges from the trajectory.

"""
from dataclasses import dataclass
from typing import Sequence, Dict, Any
from double_pendulum import DoublePendulumSimulator
from enum import Enum, auto 
from scipy.integrate import solve_ivp
import os
import numpy as np
import copy
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec

#Stablish automatically: font sizes, grid visibility, color harmony, spacing
plt.style.use("seaborn-v0_8-paper")
sns.set_theme(context="notebook", style="whitegrid", palette="viridis", font_scale=1.2)

##
#---------------- Parameters and control function -------------------
##

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

###
# --------------------- Main functions graphics --------------------------
###

def regime_summary(sol:Sequence[float], times: Sequence[float], energy: Sequence[float], position: Sequence[float], colors: Sequence[float], name: str)-> plt.figure:
    """
    Regime summary -- Position, energies and phase space
    """
    if sol is None:
        raise ValueError("Run the simulation first.")

    theta1, theta2, omega1, omega2 = sol
    time = times
    T, U, Et = energy
    x1, x2, y1, y2 = position

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
    ax2.annotate(f"Drift = {(Et[-1] - Et[0])/Et[0]:.2e}", xy=(0.05, 0.1), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round", fc="white", alpha=0.7))

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

def double_pendulum_animation(sol:Sequence[float], times:Sequence[float], energy:Sequence[float], position:Sequence[float], colors:Sequence[float], name: str)-> plt.figure:

    theta1, theta2, omega1, omega2 = sol
    
    def wrapped_theta(theta:Sequence[float])->np.ndarray:
        return (theta + np.pi) % (2*np.pi) - np.pi
    
    theta1 = wrapped_theta(theta1)
    theta2 = wrapped_theta(theta2)

    time = times
    x1, x2, y1, y2 = position
    T, U, Et = energy

    if len(time) != len(x1) or len(x1) != len(T):
        raise TypeError("Something goes wrong. The length between energy, position and sol don't fix it.")
    
    fig = plt.figure(figsize= (12, 6))
    gs = GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1], figure=fig)
    
    ax_anim = fig.add_subplot(gs[:, 0])   # Left: animation (takes both rows)
    ax_top  = fig.add_subplot(gs[0, 1])   # Right-top: PHASE SPACE
    ax_bot  = fig.add_subplot(gs[1, 1])   # Right-bottom: ANGULAR SPACE

    # --- MOTION PANEL ---

    rod1, = ax_anim.plot([], [], "k-", lw=2)
    rod2, = ax_anim.plot([], [], "k-", lw=2)
    m1, = ax_anim.plot([], [], "o", color=colors["mass1"], markersize=10)
    m2, = ax_anim.plot([], [], "o", color=colors["mass2"], markersize=10)
    ax_anim.set_title("Double pendulum motion")
    ax_anim.set_xlim(min(x2)-1, max(x2)+1)
    ax_anim.set_ylim(min(y2)-1, 1)
    ax_anim.set_xlabel("x [m]")
    ax_anim.set_ylabel("y [m]")
    ax_anim.axhline(0, color="black", lw=0.5)

    # --- PHASE SPACE PANEL ---
    ps1, = ax_top.plot([], [], "--", color=colors["mass1"], lw=2)
    ps2, = ax_top.plot([], [], "--", color=colors["mass2"], lw=2)
    dot1, = ax_top.plot([], [], "o", color=colors["mass1"])
    dot2, = ax_top.plot([], [], "o", color=colors["mass2"])

    combine_omega = omega1 + omega2
    max_y = max(combine_omega)
    min_y = min(combine_omega)

    ax_top.set_title(r"Phase space $\theta$ vs $\omega$")
    ax_top.set_xlim([-np.pi, np.pi])
    ax_top.set_ylim([min_y, max_y])
    ax_top.set_xlabel(r"$\theta$ [rad]")
    ax_top.set_ylabel(r"$\omega$ [rad/s]")

    # --- Angular plane ---
    line1, = ax_bot.plot([], [], "--", color = colors["mass1"], lw =2)
    point1, = ax_bot.plot([], [], "o", color = colors["mass1"], lw =2)
    ax_bot.set_xlabel(r"$\theta_1 [rad]$")
    ax_bot.set_ylabel(r"$\theta_2 [rad]$")
    ax_bot.set_xlim([-np.pi, np.pi])
    ax_bot.set_ylim([-np.pi, np.pi])


    # --- UPDATE FUNCTION ---
    def update(i):
        # Motion
        rod1.set_data([0, x1[i]], [0, y1[i]])
        rod2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        m1.set_data([x1[i]], [y1[i]])
        m2.set_data([x2[i]], [y2[i]])

        # Phase space
        ps1.set_data([theta1[:i]], [omega1[:i]])
        ps2.set_data([theta2[:i]], [omega2[:i]])
        dot1.set_data([theta1[i]], [omega1[i]])
        dot2.set_data([theta2[i]], [omega2[i]])

        #Phase space (thetas)
        line1.set_data([theta1[:i]], [theta2[:i]])
        point1.set_data([theta1[i]], [theta2[i]])

        return rod1, rod2, m1, m2, ps1, ps2, dot1, dot2, line1, point1
    
    plt.tight_layout()
    frame_step = 50
    anim = FuncAnimation(fig, update, frames = np.arange(0, len(time), frame_step), interval = 50, blit = False, repeat = False)
    return anim

def compute_initial_angles(params, theta1: Sequence[float], theta2: Sequence[float], directory: str = os.getcwd(), filename: str = "compute_intial_angles.npz") -> Dict[str, Any]:
    """
    Compute energy drift for:
    - varying theta1 with theta2 = 0
    - varying theta2 with theta1 = 0
    - full grid (theta1, theta2)

    Returns a dictionary with three entries:
        "theta1_scan" - "theta_1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []
        "theta2_scan" - "theta_1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []
        "grid" - "theta_1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []
    Uses:
    - Fractal
    - Errors calculus (main_error.py)

    """
    #Validation
    if len(theta1) == 0 or len(theta2) == 0:
        raise ValueError("theta1 and theta2 must each contain at least one angle")

    #Define dictionaries
    theta1_scan = {th: {"theta_1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []} for th in theta1}
    theta2_scan = {th: {"theta_1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []} for th in theta2}
    grid = {th1: {th2: {"theta_1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []} for th2 in theta2} for th1 in theta1}

    def _compute(q0: Sequence[float])->dict:
        params.theta1_0, params.theta2_0 = q0
        sim = DoublePendulumSimulator(params=params)
        results = sim.run()
        T, U, Et = results.kinetic, results.potential, results.total
        theta1, theta2, omega1, omega2 = results.y

        drift = _compute_energy_drift(Et)
        return {"theta1":theta1, "theta2": theta2, "omega1": omega1, "omega2": omega2, "kinetic": T, "potencial":U, "drift": drift}

    def _compute_energy_drift(Et: Sequence[float])->float:

        energy_drift = max(Et) - min(Et)

        return energy_drift

    # Scan theta1 (theta2 = 0)
    print("Starting with theta1_0")
    for th in theta1:
        theta1_scan[th] = _compute((th, 0))

    #Scan theta2 (theta1 = 0)
    print("Starting with theta2_0")
    for th in theta2:
        theta2_scan[th] = _compute((0.0, th))

    # Full grid (theta1, theta2)
    total = len(theta1)
    bar_len = 20

    for i, th1 in enumerate(theta1):
        progress = (i + 1) / total
        filled = int(progress * bar_len)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(f"[{bar}]  {progress*100:5.1f}%   θ₁ = {th1:.4f}", end="\r", flush=True)

        for th2 in theta2:
            grid[th1][th2] = _compute((th1, th2))

    print()  # newline after progress bar

    # -----------------------------
    # 6. Save results
    # -----------------------------
    route = os.path.join(directory, filename)
    np.savez(route, theta1_scan=theta1_scan, theta2_scan=theta2_scan, grid=grid)

    return {"theta1_scan": theta1_scan, "theta2_scan": theta2_scan, "grid": grid}

def fractal():
    return 
def numerical_jacobian(f, x, eps=1e-5):
    """Compute Jacobian numerically."""
    n = len(x)
    J = np.zeros((n, n))
    fx = f(0, x)
    for i in range(n):
        dx = np.zeros(n)
        dx[i] = eps
        J[:, i] = (f(0, x + dx) - fx) / eps
    return J

def lyapunov_spectrum(params, T=200, dt=0.01, renorm_steps=10):
    """
    Compute all 4 Lyapunov exponents of the double pendulum
    using the Benettin–Wolf algorithm.
    """

    # Independent parameter copy
    p = copy.deepcopy(params)
    sim = DoublePendulumSimulator(params=p)

    # Initial state
    x = np.array([p.theta1_0, p.theta2_0, p.omega1_0, p.omega2_0], dtype=float)

    # Initial perturbation matrix (identity)
    Q = np.eye(4)

    # Accumulated log norms
    lyap_sum = np.zeros(4)

    steps = int(T / dt)

    for step in range(steps):

        # Integrate main system
        sol = solve_ivp(sim.equations_of_motion, [0, dt], x, max_step=dt)
        x = sol.y[:, -1]

        # Compute Jacobian at new point
        J = numerical_jacobian(sim.equations_of_motion, x)

        # Evolve perturbation matrix
        Q = solve_ivp(lambda t, y: (J @ y.reshape(4,4)).flatten(),
              [0, dt], Q.flatten(),
              max_step=dt).y[:, -1].reshape(4,4)

        # Every few steps: QR decomposition
        if (step + 1) % renorm_steps == 0:
            Q, R = np.linalg.qr(Q)

            # Accumulate logs of diagonal of R
            lyap_sum += np.log(np.abs(np.diag(R)))

    # Convert sums to exponents
    total_time = steps * dt
    lyap = lyap_sum / (total_time)

    # Sort from largest to smallest
    lyap_sorted = np.sort(lyap)[::-1]

    return lyap_sorted

def lyapunov_graphics(params:Sequence[float], theta1:Sequence[float])->plt.figure:
    print("Starting the Lyapunov graphics")

    L1, L2, L3, L4 = [], [], [], []

    for i, theta in enumerate(theta1):
        bar_len = 20
        progress = (i + 1) / len(theta1)
        filled = int(progress * bar_len)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(rf"[{bar}]  {progress*100:5.1f}%   θ₁ = {theta:.4f}", end="\r", flush=True)

        params.theta1_0 = theta
        params.theta2_0 = 0.0
        params.omega1_0 = 0.0
        params.omega2_0 = 0.0

        lya = lyapunov_spectrum(params=params)

        L1.append(lya[0])
        L2.append(lya[1])
        L3.append(lya[2])
        L4.append(lya[3])

    print("\n" + "="*60)
    print("Starting the plotting")
    print("="*60)

    fig = plt.figure(figsize=(10,6))
    plt.plot(theta1, L1, label=r"$\lambda_1$ ", lw=2)
    plt.plot(theta1, L2, label=r"$\lambda_2$", lw=2)
    plt.plot(theta1, L3, label=r"$\lambda_3$", lw=2)
    plt.plot(theta1, L4, label=r"$\lambda_2$", lw=2)

    plt.axhline(0, lw=0.5, linestyle="--")
    plt.xlabel(r"$\theta_1$")
    plt.ylabel("Lyapunov exponents")
    plt.title(r"Full Lyapunov spectrum vs $\theta_1$")
    plt.legend()
    plt.xlim([0, theta1[-1]])

    return fig

def poincare():
    return 

##
# ---------------------- Position ----------------
##
def position(sol: Sequence[float],) -> np.ndarray:
    theta1 , theta2, _, _ = sol
    x1 = params.L1 * np.sin(theta1)
    x2 = x1 + params.L2 * np.sin(theta2)
    y1 = - params.L1 * np.cos(theta1)
    y2 = y1 - params.L2 * np.cos(theta2)

    return x1, x2, y1, y2
##
# ---------------------- Settle the parameters -----------------------------
##
params = Params()
sim = DoublePendulumSimulator(params = params)

# --- Graphics ---
results = sim.run()
best_method  = "Radau (implicit)"

sol = results.y
energy = results.kinetic, results.potential, results.total
times = results.t

positions = position(sol = sol)

cmap = plt.colormaps["viridis"]
color = cmap(np.linspace(0, 1, 3))
colors = {
    "mass1": cmap(0.2),
    "mass2": cmap(0.8),
    "T": "#1f77b4",
    "U": "#ff7f0e",
    "Et": "#2ca02c",
}

#fig1 =regime_summary(sol = sol, times = times, energy=energy, position=positions, colors = colors, name = "Radau (implicit)")
#fig2 = double_pendulum_animation(sol = sol, times = times, energy=energy, position=positions, colors= colors, name = "Radau (implicit)")

theta1 = np.linspace(np.deg2rad(0), np.deg2rad(45), 45)
fig = lyapunov_graphics(params=params, theta1=theta1)
plt.show()
