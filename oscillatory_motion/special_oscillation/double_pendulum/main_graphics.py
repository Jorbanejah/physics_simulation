"""
MAIN GRPAHICS: double pendulum

This code generate the following graphics such as linearizad equation as normal equation

- Regime summary (position, energies, phase space) 
- Animation
- Lyapunov coefficient


Another performance (in process):
- Poincare sections
- Trajectory + Poincaré section side-by-side: Shows how the chaotic cloud emerges from the trajectory.

"""
from dataclasses import dataclass
from typing import Sequence, Callable
from double_pendulum import DoublePendulumSimulator
from enum import Enum, auto 
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.gridspec import GridSpec

import os
import numpy as np
import copy
import matplotlib.pyplot as plt



#Stablish automatically: font sizes, grid visibility, color harmony, spacing
plt.style.use("seaborn-v0_8-paper")

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

    theta1_0: float = np.deg2rad(120)
    theta2_0: float = np.deg2rad(55)

    omega1_0: float = 0.0
    omega2_0: float = 0.0

    t_max: float = 100.0  # s
    dt: float = 0.01  # s

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
def wrapped_theta(theta:Sequence[float])->np.ndarray:
        return (theta + np.pi) % (2*np.pi) - np.pi
    
def regime_summary(sol:Sequence[float], times: Sequence[float], energy: Sequence[float], position: Sequence[float], colors: Sequence[float], name: str)-> plt.Figure:
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
    theta1 = wrapped_theta(theta1)
    theta2 = wrapped_theta(theta2)

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

def double_pendulum_animation(sol:Sequence[float], times:Sequence[float], energy:Sequence[float], position:Sequence[float], colors:Sequence[float], name: str)-> plt.Figure:

    theta1, theta2, omega1, omega2 = sol
    
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

    ax_top.set_title(r"Phase space $\theta$ vs $\omega$")
    ax_top.set_xlim([-np.pi, np.pi])
    ax_top.set_ylim(min(omega1.min(), omega2.min()), max(omega1.max(), omega2.max()))
    ax_top.set_xlabel(r"$\theta$ [rad]")
    ax_top.set_ylabel(r"$\omega$ [rad/s]")

    # --- Angular plane ---
    line1, = ax_bot.plot([], [], "--", color = colors["mass1"], lw =2)
    point1, = ax_bot.plot([], [], "o", color = colors["mass1"], lw =2)
    ax_bot.set_xlabel(r"$\theta1 [rad]$")
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
    frame_step = 5
    anim = FuncAnimation(fig, update, frames = np.arange(0, len(time), frame_step), interval = 50, blit = False, repeat = False)
    return anim

##
# --------------------- Lyapunov coefficient --------------------------
##
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

def lyapunov_graphics(params:Sequence[float], theta1:Sequence[float], color:Sequence[float], store_data: bool = False)->plt.figure:
    print("Starting the Lyapunov graphics")

    L1, L2, L3, L4 = [], [], [], []

    route = os.path.join(os.getcwd(), "Lyapunov_coeffcient.npz")
    flag = True

    for i, theta in enumerate(theta1):
        if os.path.exists(route):
            print(f"File already exists: {route}. Stopping loop.")
            break
        flag = False
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

    if store_data:
        np.savez(route, L1 = L1, L2 = L2, L3 = L3, L4 =L4)
    
    if flag:
        data = np.load(route, allow_pickle= True)
        L1 = data["L1"] 
        L2 = data["L2"] 
        L3 = data["L3"] 
        L4 = data["L4"] 

    print("\n" + "="*60)
    print("Starting the plotting")
    print("="*60)

    fig = plt.figure(figsize=(10,6))
    plt.plot(theta1, L1, label=r"$\lambda_1$ ", color = color["Ly1"], lw=2)
    plt.plot(theta1, L2, label=r"$\lambda_2$", color = color["Ly2"], lw=2)
    plt.plot(theta1, L3, label=r"$\lambda_3$", color = color["Ly3"], lw=2)
    plt.plot(theta1, L4, label=r"$\lambda_4$", color = color["Ly4"], lw=2)

    plt.axhline(0, lw=0.5, linestyle="--")
    plt.xlabel(r"$\theta1$")
    plt.ylabel("Lyapunov exponents")
    plt.title(r"Full Lyapunov spectrum vs $\theta1$")
    plt.legend()
    plt.xlim([theta1[0], theta1[-1]])

    return fig

##
# ---------------------- Position ----------------
##
def position(sol: Sequence[float], params: Sequence[float],) -> np.ndarray:
    theta1 , theta2, _, _ = sol
    x1 = params.L1 * np.sin(theta1)
    x2 = x1 + params.L2 * np.sin(theta2)
    y1 = - params.L1 * np.cos(theta1)
    y2 = y1 - params.L2 * np.cos(theta2)

    return x1, x2, y1, y2
##
# ---------------------- Store the graphics -----------------------------
##

def compute(store: bool = False)-> plt.Figure:
    """
    Compute all functions in this code

    Paramters
    ---------
    store: bool
    Control whether the figures are stored or not in the currents dictionary
    """
    params = Params()
    sim = DoublePendulumSimulator(params = params)

    # --- Graphics ---
    results = sim.run()

    sol = results.y
    energy = results.kinetic, results.potential, results.total
    times = results.t

    positions = position(sol = sol, params= params)

    cmap = plt.colormaps["viridis"]
    colors = {
        "mass1": cmap(0.2),
        "mass2": cmap(0.8),
        "T": "#1f77b4",
        "U": "#ff7f0e",
        "Et": "#2ca02c",
        "Ly1": cmap(0.3),
        "Ly2": cmap(0.4),
        "Ly3": cmap(0.5),
        "Ly4": cmap(0.6)
    }
    # Compute regme summary
    print("Starting: regime summary")
    fig1 =regime_summary(sol = sol, times = times, energy=energy, position=positions, colors = colors, name = "Radau (implicit)")
    # Animation
    print("Starting: double-pendulum animation")
    anim1 = double_pendulum_animation(sol = sol, times = times, energy=energy, position=positions, colors= colors, name = "Radau (implicit)")
    #Lyapunov coefficient
    print("Starting: Lyapunov coefficients")
    theta1 = np.linspace(np.deg2rad(-90), np.deg2rad(90), 10)
    fig3 = lyapunov_graphics(params=params, theta1=theta1, color = colors)

    if store:
        directory = os.getcwd()

        # Ensure folder exists
        fig_dir = os.path.join(directory, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        # --- Save static figures ---
        fig1.savefig(os.path.join(fig_dir, "Regime_summary.png"), dpi=300, bbox_inches='tight')

        fig3.savefig(os.path.join(fig_dir, "Lyapunov_coefficients.png"), dpi=300, bbox_inches='tight')

        # --- Save animation ---
        writer = PillowWriter(fps=30)
        anim1.save(os.path.join(fig_dir, "double_pendulum_anim.gif"), writer=writer)

    return anim1, fig1, fig3

writer = PillowWriter(fps=50)

params = Params()
sim = DoublePendulumSimulator(params = params)
results = sim.run()

sol = results.y
energy = results.kinetic, results.potential, results.total
times = results.t

positions = position(sol = sol, params= params)

cmap = plt.colormaps["viridis"]
colors = {
        "mass1": cmap(0.2),
        "mass2": cmap(0.8),
        "T": "#1f77b4",
        "U": "#ff7f0e",
        "Et": "#2ca02c",
        "Ly1": cmap(0.3),
        "Ly2": cmap(0.4),
        "Ly3": cmap(0.5),
        "Ly4": cmap(0.6)
    }

directory = os.getcwd()
fig_dir = os.path.join(directory, "figures")
anim1 = double_pendulum_animation(sol = sol, times = times, energy=energy, position=positions, colors= colors, name = "Radau (implicit)")
anim1.save(os.path.join(fig_dir, "double_pendulum_chaos_anim.gif"), writer=writer)

plt.show()