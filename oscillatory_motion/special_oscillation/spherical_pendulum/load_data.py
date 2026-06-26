import matplotlib.pyplot as plt
import numpy as np
import os
from spherical_pendulum import Spherical_Pendulum
from dataclasses import dataclass
from typing import Tuple, Sequence, Dict, Any

##
# -------------- Parameters definition -----------------------------
##
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


def compute(theta:float, phi:float, small_angle: bool)->float:

    Params.q0 = (np.deg2rad(theta), np.deg2rad(phi))

    sim = Spherical_Pendulum(params= Params, small_angle = small_angle)
    solution = sim.run()
    _, _, Et = sim.energies()

    return (np.abs(max(Et)) - np.abs(min(Et)))
    

def stored(theta_values: Sequence[float] = np.linspace(-np.pi, np.pi, 180), phi_values:Sequence[float] = np.linspace(0, np.pi, 90)) -> Dict[Any, Any, Any]:
    
    pass