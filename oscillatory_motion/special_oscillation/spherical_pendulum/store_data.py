import os
import numpy as np
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

    t:float = 20
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


def run_simulation(theta_0: float, phi_0: float, sim) ->Tuple[Any, Any, Any]:
    "This function simulate the spherical pendulum given some initial condition"

    Params.q0 = (np.deg2rad(theta_0),np.deg2rad(phi_0))

    results, runtime = sim.run(Params)
    _, _, Et = sim.energies()

    times = results["t"]
    solution = results["y"]

    return solution, times, Et

def data_structure(params : Params, theta_values:Sequence[float], phi_values: Sequence[float], small_angle:bool, method:str, directory:str = os.getcwd(), 
                   filename:str = "stored_data.npz"):
    """
    Compute energy drift for varying initial conditions.

    Parameters
    ----------
    params : Params
        Physical parameters of the pendulum.
    theta_values : Sequence[float]
        List of initial axial angles.
    phi_values : Sequence[float]
        List of initial polar angles.
    small_angles : bool
        If True, use the small-angle approximation equations.
    method: str
        You can use: RK4, RK45, Verlet (only with linearized equation) and DOP853
    directory : str
        Path to save the output file.
    filename : str
        Name of the output file.

    Returns
    -------
    Dict
        Nested dictionary containing scan results and grid results.
    """

    if len(theta_values) == 0 or len(phi_values) == 0:

        raise ValueError("theta_values or phi_values must each contain at least one angle")

    #Create directory:
    grid = {th: {phi: {} for phi in phi_values} for th in theta_values}

    #Full grid theta x phi

    print("\n" + "=" * 60)
    print(f"Building Full Grid ({method} model)")
    print("=" * 60)
    
    total_steps = len(theta_values) * len(phi_values)
    count = 0
    
    sim = Spherical_Pendulum(small_angle = small_angle, method = method)

    for th in theta_values:
        for ph in phi_values:
            count += 1
            print(f"Grid Progress: {100*count/total_steps:.1f}%", end="\r")

            grid[th][ph] = run_simulation(th, ph, sim)

    print("\nSaving data...")
    
    # Save as numpy compressed archive
    if small_angle:
        filename = "store_data_linearized.npz"
        route = os.path.join(directory, filename)
        np.savez(route, grid=grid)

    else:
        route = os.path.join(directory, filename)
        np.savez(route, grid=grid)
    

# 360 X 180
phi_values = np.arange(-180, 181, 1)

eps = 1e-8
theta_values = np.arange(-90, 91, 1)
theta_values = np.where(theta_values == 0, eps, theta_values) # Avoid the theta == 0 and therefore singularity


data_structure(params=Params, theta_values=theta_values, phi_values=phi_values,small_angle=False, method = "DOP853")