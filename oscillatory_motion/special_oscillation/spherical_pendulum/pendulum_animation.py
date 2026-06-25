
"""
En la animacion de las subplots, haz un vector unitario en la bola n(t) = r(t)/L

"""


import numpy as np
import matplotlib.pyplot as plt
from spherical_pendulum import Spherical_Pendulum
from typing import Sequence, Tuple
from matplotlib.animation import FuncAnimation, PillowWriter
from dataclasses import dataclass


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

def contour_animation(params:Params,t:Sequence[float], x:Sequence[float], y: Sequence[float], z:Sequence[float], name:str = "Projections", Flag:bool = False,):
    """
    Animated spherical pendulum with projected trajectories on the walls.
    """

    # ---- Colors ----
    cmap = plt.colormaps["viridis"]
    colors = {
        "mass": cmap(0.2),
        "xy":   cmap(0.4),
        "xz":   cmap(0.6),
        "yz":   cmap(0.8)
    }

    # ---- Figure and 3D axis ----
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection='3d')

    # ---- Projected trajectories (animated) ----
    xy_line, = ax.plot([], [], [], color=colors["xy"], lw=1.5)  # projection on XY plane
    xz_line, = ax.plot([], [], [], color=colors["xz"], lw=1.5)  # projection on XZ plane
    yz_line, = ax.plot([], [], [], color=colors["yz"], lw=1.5)  # projection on YZ plane

    # ---- Main pendulum motion ----
    rod,  = ax.plot([], [], [], "k-", lw=2)
    mass, = ax.plot([], [], [], "o", color=colors["mass"], markersize=10)

    # ---- Axis limits ----
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    z_min, z_max = min(z), max(z)

    ax.set_xlim([x_min -0.5, x_max + 0.5])
    ax.set_ylim([y_min -1, y_max + 0.5])
    ax.set_zlim([z_min -1, z_max + 1])
    # ---- Horizontal and vertical lines ------

    horizontal_x = np.linspace(x_min -0.5, x_max + 0.5, 100)
    horizontal_y = np.linspace(y_min -1, y_max + 0.5, 100)
    horizontal_z = np.linspace(z_min -1, 0, 100)

    x_line = ax.plot([horizontal_x], [0], [z_min -1], "k--", lw =0.5) 
    y_line = ax.plot([0], [horizontal_y], [z_min -1], "k--", lw =0.5) 
    z_line = ax.plot([0], [0], [horizontal_z], "k--", lw =0.5) 


    ax.set_title("Spherical Pendulum with Wall Projections")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")

    # ---- Update function ----
    def update(i):
        # XY projection (z constant)
        xy_line.set_data(x[:i], y[:i])
        xy_line.set_3d_properties(np.full(i, z_min -1))

        # XZ projection (y constant)
        xz_line.set_data(x[:i], np.full(i, y_max + 1))
        xz_line.set_3d_properties(z[:i])

        # YZ projection (x constant)
        yz_line.set_data(np.full(i, x_min -.5), y[:i])
        yz_line.set_3d_properties(z[:i])

        # Pendulum rod
        rod.set_data([0, x[i]], [0, y[i]])
        rod.set_3d_properties([0, z[i]])

        # Mass
        mass.set_data([x[i]], [y[i]])
        mass.set_3d_properties([z[i]])

        return rod, mass, xy_line, xz_line, yz_line

    # ---- Animation ----
    frame_step = 5
    anim = FuncAnimation(fig, update, frames=np.arange(0, len(t), frame_step), interval=50, blit=False, repeat=False)
    if Flag:
        import os
        directory = os.getcwd()

        # Ensure folder exists
        fig_dir = os.path.join(directory, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        # --- Save animation ---
        writer = PillowWriter(fps=30)
        anim.save(os.path.join(fig_dir, f"{name}.gif"), writer=writer)
    return anim
    return anim


def subplots_animation(sol: Sequence[float], times:Sequence[float], energy:Sequence[float], position:Sequence[float], name: str = "Subplots_animation", Flag:bool = False)-> plt.Figure:

    theta, dtheta, phi, dphi = sol
    def wrapped_theta(theta:Sequence[float])->np.ndarray:
        return (theta + np.pi) % (2*np.pi) - np.pi
    
    cmap = plt.colormaps["viridis"]
    colors = {
        "mass": cmap(0.2),
        "ps1": cmap(0.4),
        "Kinetic": cmap(0.5),
        "Potential": cmap(0.7),
        "Total": cmap(0.9)
    }

    theta = wrapped_theta(theta)
    phi = wrapped_theta(phi)

    time = times
    x, y, z = position
    T, U, Et = energy

    if len(time) != len(x) or len(x) != len(T):
        raise TypeError("Something goes wrong. The length between energy, position and sol don't fix it.")
    

    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize= (12, 6))
    gs = GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1], figure=fig)

    from mpl_toolkits.mplot3d import Axes3D

    ax_anim = fig.add_subplot(gs[:, 0], projection="3d")
    ax_top  = fig.add_subplot(gs[0, 1])   # Right-top: PHASE SPACE
    ax_bot  = fig.add_subplot(gs[1, 1])   # Right-bottom: ANGULAR SPACE

    # --- MOTION PANEL ---

    rod1, = ax_anim.plot([], [], [], "k-", lw=2)
    m, = ax_anim.plot([], [], [], "o", color=colors["mass"], markersize=10)

    ax_anim.set_title("Spherical motion")
    ax_anim.set_xlim(min(x)-0.5, max(x)+0.5)
    ax_anim.set_ylim(min(y)-1, max(y)+0.5)
    ax_anim.set_zlim(min(z)-1, max(z)+1)

    ax_anim.set_xlabel("x [m]")
    ax_anim.set_ylabel("y [m]")
    ax_anim.set_zlabel("z [m]")

    # Horizontal lines
    horizontal_x = np.linspace(min(x) -0.5, max(x) + 0.5, 100)
    horizontal_y = np.linspace(min(y) -1, max(y) + 0.5, 100)
    horizontal_z = np.linspace(min(z) -1, 0, 100)

    x_line = ax_anim.plot([horizontal_x], [0], [min(z) -1], "k--", lw =0.5) 
    y_line = ax_anim.plot([0], [horizontal_y], [min(z) -1], "k--", lw =0.5) 
    z_line = ax_anim.plot([0], [0], [horizontal_z], "k--", lw =0.5)

    # --- Angular plane ---
    ps1, = ax_top.plot([], [], "--", color=colors["ps1"], lw=2)
    dot1, = ax_top.plot([], [], "o", color=colors["mass"])

    ax_top.set_title(r"Phase space $\theta$ vs $\phi$")
    ax_top.set_xlim([-np.pi, np.pi])
    ax_top.set_ylim([-np.pi, np.pi])
    ax_top.set_xlabel(r"$\theta$ [rad]")
    ax_top.set_ylabel(r"$\phi$ [rad]")

    # --- Angular plane ---
    x_values = np.arange(0, len(time))
    line1, = ax_bot.plot([], [], "-", color = colors["Kinetic"], lw =2, label= "Kinetic")
    line2, = ax_bot.plot([], [], "-", color = colors["Potential"], lw =2, label = "Potential")
    line3, = ax_bot.plot([], [], "-", color = colors["Total"], lw =2, label = "Total")

    ax_bot.set_title(r"Energies")
    ax_bot.set_xlabel(r"Time $x$ Time step")
    ax_bot.set_ylabel(r"Energy [E unit]")
    ax_bot.set_xlim([0, len(time)])
    ax_bot.set_ylim([min(U) -1, max(T) +1 ])
    ax_bot.legend(loc = "center right")


    # --- UPDATE FUNCTION ---
    def update(i):
        # Motion
        rod1.set_data([0, x[i]], [0, y[i]])
        rod1.set_3d_properties([0, z[i]])

        m.set_data([x[i]], [y[i]])
        m.set_3d_properties([z[i]])

        # Phase space (angular coordenates)
        ps1.set_data([theta[:i]], [phi[:i]])
        dot1.set_data([theta[i]], [phi[i]])

        #Energies
        line1.set_data([x_values[:i]], [T[:i]])
        line2.set_data([x_values[:i]], [U[:i]])
        line3.set_data([x_values[:i]], [Et[:i]])

        return rod1, m, ps1, dot1, line1, line2, line3,
    
    plt.tight_layout()
    frame_step = 5
    anim = FuncAnimation(fig, update, frames = np.arange(0, len(time), frame_step), interval = 50, blit = False, repeat = False)
    
    if Flag:
        import os
        directory = os.getcwd()

        # Ensure folder exists
        fig_dir = os.path.join(directory, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        # --- Save animation ---
        writer = PillowWriter(fps=30)
        anim.save(os.path.join(fig_dir, f"{name}.gif"), writer=writer)
    return anim


params = Params()

sim = Spherical_Pendulum(Params, small_angle= False, method="Rk4")

results = sim.run()

solution = results["y"]
times = results["t"]

x, y, z = sim.transform()
energies = sim.energies()

#anim1 = contour_animation(params=Params, t = times, x= x, y= y, z =z, Flag=True)
anim2 = subplots_animation(sol = solution, times = times, energy=energies, position=(x, y, z), Flag = True)
