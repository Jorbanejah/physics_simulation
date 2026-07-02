"""Spherical pendulum simulation.

The theta angle is measured from the vertical axis (polar angle), while the phi angle is measured from the horizontal axis
(azimuthal angle). The full model is integrated from the Lagrange equations.

Description:
Consider a compact mass (m) on the end of an inextansable rod of length L. This mass is free to move in any direction.
We can define our coordinate system as the fixed point at the end of the rod. Then, through cartesian equation we can define the spherical 
coordinates.

The following code calculate the trajectory with four different numerical method: {RK4 , RK45, DOP, Verlet}. Such as the spherical pendulum equation or approximation the code output it will be:

-> run() -> a sol solution with times and angular coodernate
-> transform() -> a Tuple with cartesian coordenates
-> energies() -> a Tuple with kinetic, potential and total energy.


Code by: Jorge Orbaneja Huerta
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Callable

import time
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp


#Define an alias where we will store the current state: [theta, dtheta, phi, dphi]
State = np.ndarray

Dynamics = Callable[[float, State, "Params"], State] # Input (Time, State, Params (not define yet)) ----> Output (State)


###
# -------------Params and functions that controls the physics --------------
###


# Control if the current initial state have two parameters q0 = [theta0, phi0] and dq0 = [dtheta, dphi]
def _as_pair(values: Sequence[float], name: str) -> tuple[float, float]:
    
    if len(values) != 2:
        raise ValueError(f"{name} must contain exactly two values")
    return (float(values[0]), float(values[1]))


def _time_grid(duration: float, dt: float) -> np.ndarray:
    """Return integration times that always include the final time"""

    if duration <= 0.0:
        raise ValueError("Simulation time must be positive.")
    if dt <= 0.0:
        raise ValueError("Time step must be positive.")
    
    times = np.arange(0.0, duration, dt, dtype = float)
    if times.size == 0 or not np.isclose(times[-1], duration):
        times = np.append(times, duration)
    else:
        times[-1] = duration

    return times

@dataclass(slots=True)
class Params:
    "Phisical parameters and initial condition for the pendulum"
    g: float = 9.81

    m: float = 1.0
    L:float = 2.0

    q0: tuple[float, float]= (np.deg2rad(10.0), np.deg2rad(0.0))
    dq0: tuple[float, float] = (0.0, 0.0)

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

###
# ------------------------- Integrators --------------------------
###

def velocity_verlet(f: Dynamics, t: float, dt: float, y: State, params: Params) -> State:
    """Velocity Verlet step for systems whose acceleration depends on position."""
    y = np.asarray(y, dtype=float) #Convert the input to an array
    angle = y[:2] # theta, phi
    d_angle = y[2:] # dtheta, dphi

    acceleration = f(t, y, params)[2:] # Returns only dtheta, dphi
    angle_new = angle + d_angle * dt + 0.5 * acceleration * dt**2

    # In the small-angle model acceleration is position-only, so the old
    # velocity is sufficient for evaluating the next acceleration.

    predicted_state = np.array([angle_new[0], angle_new[1], d_angle[0], d_angle[1]])

    acceleration_new = f(t + dt, predicted_state, params)[2:]

    d_angle_new = d_angle + 0.5 * (acceleration + acceleration_new) * dt

    return np.array([angle_new[0], angle_new[1], d_angle_new[0], d_angle_new[1]])


def rk4(f: Dynamics, t: float, dt: float, y: State, p: Params) -> State:
    """One fourth-order Runge-Kutta step."""
    y = np.asarray(y, dtype=float)
    k1 = f(t, y, p)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1, p)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2, p)
    k4 = f(t + dt, y + dt * k3, p)
    return y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

###
# ----------------------Equation------------------------
###

def spherical_pendulum_equation(t: float, y:State, p:Params) -> State:
    """Full nonlinear differential spherical-pendulum equations

    State vector: [theta, phi, dtheta, dphi].

    """

    del t

    L, g = p.L, p.g

    theta, phi, dtheta, dphi = np.asarray(y, dtype = float)

    tan_theta = np.tan(theta)

    # avoid singularity at theta = ±pi/2
    if abs(tan_theta) < 1e-8:
        tan_theta = np.sign(tan_theta) * 1e-8

    d2phi = -2 * dtheta * dphi / tan_theta

    d2theta = (dphi**2 * np.sin(2 * theta))/ 2 - g/L * np.sin(theta)


    return np.array([dtheta, dphi, d2theta, d2phi])

def spherical_pendulum_approx(t:float, y: State, p:Params)-> State:
    """Full linear differential spherical-pendulum equation
    
    state vector = [theta, phi, dtheta, dphi]
    """
    del t

    L, g = p.L, p.g
    
    theta, phi, dtheta, dphi = np.asarray(y, dtype=float)

    d2phi = -2 * dtheta * dphi / theta
    d2theta = (dphi**2 * 2 * theta) / 2 - g/L * theta

    return np.array([dtheta, dphi, d2theta, d2phi])

###
# --------------- Main Class -----------------
###

class Spherical_Pendulum:
    """Numerical model of a Spherical pendulum."""

    _SOLVE_IVP_METHODS = {"RK45", "DOP853"}

    def __init__(self, small_angle: bool = False, method: str = "RK4",) -> None:
        '''
        Parameters
        -----------

        Params: Dict
            spherical initial conditions and characteristhic

        Small_angle: bool  
            default FALSE

        Method: string
            numerical solver. Default: RK4
        '''
        
        self.method =  method
        self.small_angle = small_angle

        self.sol = None
        self.t =  None
        self.y =  None

    def _dynamics(self) -> Dynamics:
        return spherical_pendulum_approx if self.small_angle else spherical_pendulum_equation
    
    def _initial_state(self) -> State:
        return np.array([self.params.q0[0], self.params.q0[1], self.params.dq0[0], self.params.dq0[1],], dtype = float,)
    
    def run(self, params:Params):
        """Run the simulation using the configured integration method"""

        self.params = params

        self.sol = None
        self.y = None
        self.t = None

        method = self.method.strip().upper()
        dynamics = self._dynamics()
        y0 = self._initial_state()
        t_eval = _time_grid(self.params.t, self.params.dt)

        if method == "RK4":

            y_history = np.empty((4, t_eval.size), dtype = float)
            y_history[:, 0] = y0
            current_y = y0

            start = time.perf_counter()
            for index, (t_start, t_end) in enumerate(zip(t_eval[:-1], t_eval[1:]), start=1):

                current_y = rk4(dynamics, t_start, t_end - t_start, current_y, self.params,)
                y_history[:, index] = current_y

            runtime =  time.perf_counter() - start

            self.sol = {"t": t_eval, "y": y_history, "method": method}
            self.t = t_eval
            self.y = y_history

        elif method == "Verlet":

            if not self.small_angle:
                raise ValueError("This method only will work properly with small angle." \
                "Use RK45 OR RK4 for the full nonlinear spherical pendulum")

            y_history = np.empty((4, t_eval.size), dtype = float)
            y_history[:, 0] = y0
            current_y = y0

            start =  time.perf_counter()

            for index, (t_start, t_end) in enumerate(zip(t_eval[:-1], t_eval[1:]), start=1):

                current_y = velocity_verlet(dynamics, t_start, t_end - t_start, current_y, self.params,)
                y_history[:, index] = current_y
                
            runtime =  time.perf_counter() - start

            self.sol = {"t": t_eval, "y": y_history, "method": method}
            self.t = t_eval
            self.y = y_history

        elif method in self._SOLVE_IVP_METHODS:

            start =  time.perf_counter()
            sol = solve_ivp(dynamics, (0.0, self.params.t), y0, args=(self.params,), t_eval=t_eval, method=method, rtol=1e-9, atol=1e-11,)

            runtime =  time.perf_counter() - start

            if not sol.success:
                raise RuntimeError(f"Integration failed: {sol.message}")

            self.sol = sol
            self.t = sol.t
            self.y = sol.y

        else:
            valid = ", ".join(sorted([*self._SOLVE_IVP_METHODS, "RK4", "Verlet"]))
            raise ValueError(f"Unknown integration method {self.method!r}. Use one of: {valid}.")

        return self.sol,  runtime

    def transform(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Convert angular coordenates to Cartesian coordinates"""
        if self.y is None:
            raise ValueError("Run the simulation first.")
        
        theta, phi = self.y[0], self.y[1]
        L = self.params.L

        x = L * np.sin(theta) * np.cos(phi)
        y = L * np.sin(theta) * np.sin(phi)
        z = - L * np.cos(theta)

        return x, y, z
    
    def energies(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:

        if self.y is None:
            raise ValueError("Run the simulation first")
        
        theta, _, dtheta, dphi = self.y[0], self.y[1], self.y[2], self.y[3]

        m, L, g = self.params.m, self.params.L, self.params.g

        #Kinetic and potencial energy
        T = 0.5 * m * L **2 * (dtheta**2  + dphi**2 * np.sin(theta) **2)
        U = - m * g * L * np.cos(theta)
        

        return T, U,  T + U