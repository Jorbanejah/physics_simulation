"""
Spherical pendulum - main graphics:

Main graphics will be run by both: linearized and normal simulation.

- Time vs omega, Time vs phi - Done
- Phase space (4x4) mix (theta, phi, dtheta, dphi) - Done
- Varying initial condition:
    - Phase portrait: nutation, precession
    - Resonance ratio vs dtheta/dphi

Numerical graphics:
- Comparison method.
- Energy drift colormap.
- Runtime.
- Convergence/stability.

Special grahics:
- Invariant-torus reconstruction. Use delay embedding or Fourier decomposition to visualize the torus in 3D.
- Frequency-map analysis (Laskar). Plot frequency drift to detect weak chaos. (MIRA EN FAVORITOS DE GOOGLE)
- Action-angle coordinate plots. If you compute approximate actions, J_theta,J_phi, plot trajectories in action space.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from spherical_pendulum import Spherical_Pendulum
from dataclasses import dataclass
from typing import Tuple, Sequence, Dict, Any

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

def compute(params:Params, method:str = "Rk4", small_angle: bool = False)-> Tuple[float, float, float]:
    "Compute the spherical pendulum simulation"

    params = Params()
    sim = Spherical_Pendulum(params= Params, small_angle=small_angle, method= method)

    solution = sim.run()
    energies = sim.energies()
    cartesian = sim.transform()

    return solution, energies, cartesian

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
# ------------------------- Main Graphics --------------------------
##

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
        A figure containing all pairwise phase‑space subplots.

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

    fig, axes = plt.subplots(nrows= N, ncols= N-1, figsize = (12, 8), tight_layout =True)
    ax = axes.ravel() #Faltten - 1D plots

    # Create a symbol map to transform each letter into a symbol-Greek letter:
    
    Greek_letter = {
        "theta": r"$\theta$",
        "phi": r"$\phi$",
        "dtheta": r"$\dot{\theta}$",
        "dphi": r"$\dot{\phi}$"
    }


    k = 0

    for var in variable_list:

        try: 
            variable_1 = params[var]

        except:
            raise ValueError(f"Missing parameter for variable {var}")
        
        for vari in variable_list:

            if var == vari:
                continue
            
            try:
                variable_2 = params[vari]

            except:
                raise TypeError(f"Missing parameter for variable {vari}")
            
            ax[k].plot(variable_1, variable_2, "--", lw = 2)
            
            try:
                ax[k].set_xlabel(rf"{Greek_letter[var]}")
                ax[k].set_ylabel(rf"{Greek_letter[vari]}")

            except:
                raise ValueError(f"Variable {var} or {vari} has no Greek symbol mapping")
            k +=1

    return fig
            
            
def nutation_precession(times:Sequence[float], dtheta:Sequence[float], dphi:Sequence[float]) -> plt.Figure:

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


solution, energies, cartesian = compute(params=Params)

theta, dtheta, phi, dphi = solution["y"]
times = solution["t"]
fig = phase_space(variable_list=["theta", "phi", "dtheta", "dphi"], theta = theta, phi = phi, dtheta = dtheta, dphi = dphi)
#fig2 = nutation_precession(times = times, dtheta = dtheta, dphi= dphi)

plt.show()



