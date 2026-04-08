from dataclasses import dataclass
import numpy as np
from Driven_oscillation import DrivenOscillation
import manim as mn

@dataclass
class DrivenOscillationParameter:
    #Initial condition
    q0: float = np.deg2rad(30)
    dq0: float = 0

    #Innate parameter
    L: float = 1
    m: float = 2
    gamma: float = 0.4

    #External force
    F0: float = 2
    F0_external: str = 'cos'
    omega: float = np.deg2rad(30)

    #Time
    t: int = 20
    dt: float = 0.01

    #System
    system: str = 'nonlinear'


class DynamicalSystem:
    """
    Encapsulates numerical simulation for a sweep of the forcing amplitude.
    """

    def __init__(self, params=None, alphas=None, methods=None):
        self.params = params or DrivenOscillationParameter()
        self.alphas = alphas or np.linspace(0.1, 3, 50)
        self.methods = methods or ["rk4", "CN", "Verlet"]

        # Allocate storage
        self.phase_portrait = {
            method: {"q": {}, "dq": {}} for method in self.methods
        }

    def run_parameter_sweep(self):
        """
        Compute phase portraits for all alpha values.
        """
        for alpha in self.alphas:
            # Update F0 = alpha * m
            F0 = alpha * self.params.m

            osc = DrivenOscillation(
                q0=self.params.q0,
                dq0=self.params.dq0,
                m=self.params.m,
                gamma=self.params.gamma,
                F0=F0,
                omega=self.params.omega,
                system=self.params.system,
                L=self.params.L,
            )

            model = osc.run()

            for method in self.methods:
                q = model.history[method]["q"]
                dq = model.history[method]["dq"]

                self.phase_portrait[method]["q"][alpha] = q
                self.phase_portrait[method]["dq"][alpha] = dq

        return self.phase_portrait


class Scene():