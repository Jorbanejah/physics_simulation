"""
Spherical pendulum - main graphics:

Main graphics can be run by both: linearized and normal simulation.

- Time vs omega, Time vs phi - Done
- Phase space (4x4) mix (theta, phi, dtheta, dphi) - Done
- Varying initial condition:
    - Phase portrait: nutation, precession

"""

import matplotlib.pyplot as plt
import numpy as np
import os
from spherical_pendulum import Spherical_Pendulum
from dataclasses import dataclass
from typing import Tuple, Sequence, Dict

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

def run_simulation(params:Params, method:str = "Rk4", small_angle: bool = False)-> Tuple[float, float, float]:
    "Compute the spherical pendulum simulation"

    params = Params()
    sim = Spherical_Pendulum(small_angle=small_angle, method= method)

    solution, _ = sim.run(params= Params)
    energies = sim.energies()
    cartesian = sim.transform()

    return solution, energies, cartesian

def upload_data(name:str = "loaded_data.npz", directory: str = os.getcwd()) -> Dict:

    path = os.path.join(directory, name)

    if os.path.exists(path):
        
        data = np.load(path, allow_pickle= True)

        grid =data["grid"].item()

        return grid

    else:
        raise FileNotFoundError("Compute first the load_data.py file")

##
# ------------------------- Main Graphics --------------------------
##

def wrapped_theta(angle:Sequence[float])->np.ndarray:
        return (angle + np.pi) % (2*np.pi) - np.pi

def phase_space(variable_list: list[str], **kwargs)->plt.Figure:

    """
    Plot all pairwise combinations of the given variables (phase-space style).

    Description
    -----------
    Given a list of variable names and corresponding data arrays passed as 
    keyword arguments, this function generates a grid of subplots containing 
    every possible ordered pair (x, y) of the variables, excluding pairs where 
    a variable is plotted against itself.

    For N variables, the output is a (N-1) x (N-1) grid containing N(N-1) 
    phase-space plots.

    Parameters
    ----------
    variable_list : list[str]
        List of variable names. Each name must match a keyword argument.
    
    **kwargs :
        Data arrays corresponding to each variable name in `variable_list`.

    Returns
    -------
    matplotlib.figure.Figure
        A figure containing all pairwise phase-space subplots.

    Examples
    --------
    >>> phase_space(
    ...     ["theta", "phi", "dtheta"],
    ...     theta=theta,
    ...     phi=phi,
    ...     dtheta=dtheta
    ... )
    This produces a 2x2 grid with all ordered combinations except self-pairs.
    """

    if len(variable_list) <= 1:
        raise ValueError("The 'variable' parameter must be, at least, length two")
    
    params = kwargs
    if len(variable_list) != len(params):
        raise ValueError("The parameters number given must be the same as variable")
    
    N = len(variable_list) # THe total subplots will be N x (N-1)

    import math
    fig, axes = plt.subplots(nrows= int(math.trunc((N * N - N)/4) -1), ncols= int(math.trunc((N * N - N)/4) + 1), figsize = (12, 8), tight_layout =True)


    ax = axes.ravel() #Faltten - 1D plots

    # Create a symbol map to transform each letter into a symbol-Greek letter:
    
    Greek_letter = {
        "theta": r"$\theta$",
        "phi": r"$\phi$",
        "dtheta": r"$\dot{\theta}$",
        "dphi": r"$\dot{\phi}$"
    }

    k = 0
    stored =[] #This varaible is defined to stored the var parameter 
    for var in variable_list:

        #Store the parameter
        stored.append(var)
        try: 
            variable_1 = params[var]

        except:
            raise ValueError(f"Missing parameter for variable {var}")
        
        for vari in variable_list:
            
            if vari in stored:
                continue
            
            try:
                variable_2 = params[vari]

            except:
                raise TypeError(f"Missing parameter for variable {vari}")
            
            if var == "theta" or var == "phi":
                variable_1 = wrapped_theta(angle = variable_1)
            if vari == "theta" or vari == "phi":
                variable_2 = wrapped_theta(angle = variable_2)

            ax[k].plot(variable_1, variable_2, "--", lw = 2)
            
            try:
                ax[k].set_xlabel(rf"{Greek_letter[var]}")
                ax[k].set_ylabel(rf"{Greek_letter[vari]}")

            except:
                raise ValueError(f"Variable {var} or {vari} has no Greek symbol mapping")
            
            k +=1

    return fig
            
            
def angular_velocity_time(times:Sequence[float], dtheta:Sequence[float], dphi:Sequence[float]) -> plt.Figure:

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols =2, tight_layout = True, figsize = (10,6))

    if len(times) != len(dtheta) or len(times) != len(dphi):
        raise TypeError("Something goes wrong. Please run again the simulation. If the problem persist, fix it :)")
    
    cmap= plt.colormaps["viridis"]
    ax1.plot(times, dtheta, color = cmap(0.2))
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel(r"$\dot{\theta} [rad/s]$")
    ax1.set_title(r"$\dot{\theta}$ vs Time")

    ax2.plot(times, dphi, color = cmap(0.5))
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel(r"$\dot{\phi} [rad/s]$")
    ax2.set_title(r"$\dot{\phi}$ vs Time")

    return fig

def nutation_precession(grid:Sequence[float], number: int) -> Tuple[plt.Figure, plt.Figure]:

    """
    The function shows how the spherical pendulum describe different phase space (theta vs phi) through random initial condition.

    Parameter:
    -------------
    grid: Sequence[float]
        store data from stored_data.py
    number: int
        it controls how many random initial condition will be taken
    """
    #Extract the values and the inex
    theta_vals = sorted(grid.keys())
    phi_vals = sorted(grid[theta_vals[0]].keys())

    thetas = np.random.choice(theta_vals, size=number)
    phis = np.random.choice(phi_vals, size = number)

    #Extract the current trajectory for each index

    theta = []
    phi = []
    dtheta =[]
    dphi = []
 
    for i, (th, ph) in enumerate(zip(thetas, phis)):

        solution, _, _ = grid[th][ph]
      
        theta.append(wrapped_theta(solution[:, 0]))   # theta(t)
        phi.append(wrapped_theta(solution[:, 1]))     # phi(t)

        dtheta.append(solution[:, 2])                 # dtheta/dt
        dphi.append(solution[:, 3])                   # dphi/dt


    fig = plt.figure(figsize = (10,6))

    ax = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)

    cmap = plt.colormaps["viridis"]

    for i in range(number):
        color = i/number
        ax.plot(theta[i], dtheta[i], color = cmap(color), linewidth = 1, label = rf"$\theta_0 ={thetas[i]}, \phi_0 = {phis[i]}$")
        ax1.plot(phi[i], dphi[i], color = cmap(color), linewidth = 1, label = rf"$\theta_0 ={thetas[i]}, \phi_0 = {phis[i]}$")

    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\phi$")
    ax.set_xlim(-np.pi/2, np.pi/2)
    ax.set_ylim(-np.pi/2, np.pi/2)

    ax1.set_xlim(-np.pi/2, np.pi/2)
    ax1.set_ylim(-np.pi/2, np.pi/2)
    ax1.set_xlabel(r"$\dot{\theta}$")
    ax1.set_ylabel(r"$\dot{\phi}$")
    ax.legend()
    ax1.legend()
    fig.suptitle("Phase space with continuos colormap")
    plt.tight_layout()

    return fig

##
# -----------------------------COMPUTE GRAPHICS --------------------------------
##

def compute(flag:bool = False, stored: bool = False):
    """
    The function is desgined to run the whole error_graphics.py file. The current parameters are:

    Parameters
    -------------

    flag: bool
        the parameter controls whether the simulation will be linearized or nor. Default = False
    stored: bool
        the paremeter controls whether the current grpahics are stored or not. Default = False
    
    Friendly reminder
    -----------------
    You can change all parameters inside Parameters class: Time, length, mass...
    """
    params = Params
    results = run_simulation(params= params, method = "DOP853", small_angle=flag)

    solution, _, _ = results
    sol = solution["y"]
    times = solution["t"]
    theta = sol[:, 0]
    dtheta =sol[:, 1]
    phi = sol[:, 2]
    dphi = sol[:, 3]

    print("Starting with figures:")
    
    phase = phase_space(variable_list=["theta", "phi", "dtheta", "dphi"],theta =theta, phi = phi, dtheta = dtheta, dphi = dphi)
    angular_valocity = angular_velocity_time(times = times, dtheta= dtheta, dphi=dphi)
    
    print("Starting data files and stored the figures")

    directory = os.getcwd()
    route = os.path.join(directory, "figures")
    os.makedirs(route, exist_ok=True)

    if flag:
            
        grid = upload_data(name = "linearized_stored_data.npz")
            
        nut_press = nutation_precession(grid = grid, numeber = 4)

        if stored:
            phase.savefig(fname= os.path.join(route, "phase_space_linearized.png"), dpi = 300, bbox_inches = "tight")
            angular_valocity.savefig(fname= os.path.join(route, "angular_velocity_linearized.png"), dpi = 300, bbox_inches = "tight")
            nut_press.savefig(fname= os.path.join(route, "nutation_preccesion_linearized.png"), dpi = 300, bbox_inches = "tight")

        else:
            plt.show()

    else:

        grid = upload_data(name = "stored_data.npz")

        nut_press = nutation_precession(grid = grid, numeber = 4)
            
        if stored:

            phase.savefig(fname= os.path.join(route, "phase_space.png"), dpi = 300, bbox_inches = "tight")
            angular_valocity.savefig(fname= os.path.join(route, "angular_velocity.png"), dpi = 300, bbox_inches = "tight")
            nut_press.savefig(fname= os.path.join(route, "nutation_preccesion.png"), dpi = 300, bbox_inches = "tight")
        else:
            plt.show()

if __name__ == "__main__":
    grid =upload_data(name = "stored_data.npz")
    nutation_precession(grid, number =4)
    #compute(flag = False, stored = True)

plt.show()



