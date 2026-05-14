"""Double pendulum simulation.

The angles are measured from the vertical axis. The full model is integrated
from the Lagrange equations using the physical mass matrix, which keeps the
signs and coupling terms explicit.


Description:

"""

from __future__ import annotations # It lets me define varaible types that it has not already defined. This kind of function goes with "_"

from collections.abc import Sequence # This is a command abstract type would let me represent whatever structure (list, tuple, string...) without involve it.
from dataclasses import dataclass #With that, python will generate __init__/ __repr__/ __eq__
from typing import Callable # It lets me describe the function type. It is use for those functions which pass as arguments

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


State = np.ndarray #Define an alias where we will store the current state: [th1, w1, th2, w2]
Dynamics = Callable[[float, State, "Params"], State] #It is a function whose inputs are: (time (float), state, Params), and its ouput is State


###
#------- Parameters and functions where control the physics of the system -------
#
# _as_pair() -----> The function control whether the current state has two values for q and dq
#_time_grid() -----> The function control whether both the time and the time step is positive
# Class() ----> where you define the main paramenter showing whether or not are positive
#
###

def _as_pair(values: Sequence[float], name: str) -> tuple[float, float]:
    if len(values) != 2:
        raise ValueError(f"{name} must contain exactly two values.")
    return float(values[0]), float(values[1])


def _time_grid(duration: float, dt: float) -> np.ndarray:
    """Return integration times that always include the final time."""
    if duration <= 0.0:
        raise ValueError("Simulation time must be positive.")
    if dt <= 0.0:
        raise ValueError("Time step must be positive.")

    times = np.arange(0.0, duration, dt, dtype=float)
    
    #Force exact final time
    times[-1] = duration
    return times


@dataclass(slots=True) #Attributes are stored in a compact C‑level array, however, you cannot describe new variables
class Params:
    """Physical parameters and initial conditions for the pendulum."""

    g: float = 9.81  # m/s^2

    m1: float = 1.0  # kg
    m2: float = 1.5  # kg
    L1: float = 1.0  # m
    L2: float = 2.0  # m

    q0: tuple[float, float] = (np.deg2rad(10.0), np.deg2rad(10.0))  # rad
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


###
# -------------- Integrators ----------------
###

def velocity_verlet(f: Dynamics, t: float, y: State, dt: float, params: Params) -> State:
    """Velocity Verlet step for systems whose acceleration depends on position."""
    y = np.asarray(y, dtype=float) #Convert the input to an ndarray
    theta = y[:2] # th1, th2
    omega = y[2:] # w1, w2

    acceleration = f(t, y, params)[2:] # Returns only w1, w2
    theta_new = theta + omega * dt + 0.5 * acceleration * dt**2

    # In the small-angle model acceleration is position-only, so the old
    # velocity is sufficient for evaluating the next acceleration.

    predicted_state = np.array([theta_new[0], theta_new[1], omega[0], omega[1]])

    acceleration_new = f(t + dt, predicted_state, params)[2:]

    omega_new = omega + 0.5 * (acceleration + acceleration_new) * dt

    return np.array([theta_new[0], theta_new[1], omega_new[0], omega_new[1]])


def rk4(f: Dynamics, t: float, dt: float, y: State, p: Params) -> State:
    """One fourth-order Runge-Kutta step."""
    y = np.asarray(y, dtype=float)
    k1 = f(t, y, p)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1, p)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2, p)
    k4 = f(t + dt, y + dt * k3, p)
    return y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


###
#-----------------Equations-----------------
###

def equation_double_pendulum(t: float, y: State, p: Params) -> State:
    """Full nonlinear double-pendulum state derivative.

    State vector: [theta_1, theta_2, omega_1, omega_2].

    """
    del t
    theta1, theta2, omega1, omega2 = np.asarray(y, dtype=float)
    delta = theta1 - theta2
    sin_delta = np.sin(delta)
    cos_delta = np.cos(delta)

    m1, m2 = p.m1, p.m2
    L1, L2 = p.L1, p.L2
    g = p.g

    mass_matrix = np.array(
        [
            [(m1 + m2) * L1**2, m2 * L1 * L2 * cos_delta],
            [m2 * L1 * L2 * cos_delta, m2 * L2**2],
        ],
        dtype=float,
    )
    forcing = np.array(
        [
            -(m1 + m2) * g * L1 * np.sin(theta1) - m2 * L1 * L2 * omega2**2 * sin_delta,
            m2 * L1 * L2 * omega1**2 * sin_delta - m2 * g * L2 * np.sin(theta2),
        ],
        dtype=float,
    )
    # Matrix: 
    #      (A_11, A_12) (alpha1) = (b1)
    #      (A_21, A_22) (alpha2) = (b2)
    # A X = B ----> X = A^-1 * B

    alpha1, alpha2 = np.linalg.solve(mass_matrix, forcing)
    return np.array([omega1, omega2, alpha1, alpha2])


def equation_approx(t: float, y: State, p: Params) -> State:
    """Small-angle linearized double-pendulum state derivative.
    |theta_1|, |theta_2| << 1
    """
    del t
    theta1, theta2, omega1, omega2 = np.asarray(y, dtype=float)

    m1, m2 = p.m1, p.m2
    L1, L2 = p.L1, p.L2
    g = p.g

    alpha1 = -g * ((m1 + m2) * theta1 - m2 * theta2) / (m1 * L1)
    alpha2 = g * (m1 + m2) * (theta1 - theta2) / (m1 * L2)

    return np.array([omega1, omega2, alpha1, alpha2])

###
# ----------- Main Class -----------
###

class DoublePendulum:
    """Numerical model of a planar double pendulum."""

    _SOLVE_IVP_METHODS = {"RK45", "DOP853"}

    def __init__(self, params: Params, small_angle: bool = False, method: str = "RK45",) -> None:
        self.params = params
        self.small_angle = small_angle
        self.method = method
        self.sol = None
        self.t: np.ndarray | None = None # The varaible can be either a numpy ndarray or None
        self.y: np.ndarray | None = None

    def _dynamics(self) -> Dynamics:
        return equation_approx if self.small_angle else equation_double_pendulum

    def _initial_state(self) -> State:
        return np.array([self.params.q0[0], self.params.q0[1], self.params.dq0[0], self.params.dq0[1],], dtype=float,)

    def run(self):
        """Run the simulation using the configured integration method."""
        self.sol = None
        self.t = None
        self.y = None

        method = self.method.strip().upper()
        dynamics = self._dynamics()
        y0 = self._initial_state()
        t_eval = _time_grid(self.params.t, self.params.dt)

        if method == "RK4":

            y_history = np.empty((4, t_eval.size), dtype=float)
            y_history[:, 0] = y0
            current_y = y0

            for index, (t_start, t_end) in enumerate(zip(t_eval[:-1], t_eval[1:]), start=1):

                current_y = rk4(dynamics, t_start, t_end - t_start, current_y, self.params)
                y_history[:, index] = current_y

            self.sol = {"t": t_eval, "y": y_history, "method": method}
            self.t = t_eval
            self.y = y_history

        elif method == "VERLET":
            if not self.small_angle:

                raise ValueError(
                    "Velocity Verlet is only valid here for the small-angle model. "
                    "Use RK45 or RK4 for the full nonlinear pendulum."
                )

            y_history = np.empty((4, t_eval.size), dtype=float)
            y_history[:, 0] = y0
            current_y = y0

            for index, (t_start, t_end) in enumerate(zip(t_eval[:-1], t_eval[1:]), start=1):

                current_y = velocity_verlet(dynamics, t_start, current_y, t_end - t_start, self.params,)
                y_history[:, index] = current_y

            self.sol = {"t": t_eval, "y": y_history, "method": method}
            self.t = t_eval
            self.y = y_history

        elif method in self._SOLVE_IVP_METHODS:
            sol = solve_ivp(dynamics, (0.0, self.params.t), y0, args=(self.params,), t_eval=t_eval, method=method, rtol=1e-9, atol=1e-11,)

            if not sol.success:
                raise RuntimeError(f"Integration failed: {sol.message}")

            self.sol = sol
            self.t = sol.t
            self.y = sol.y

        else:
            valid = ", ".join(sorted([*self._SOLVE_IVP_METHODS, "RK4", "Verlet"]))
            raise ValueError(f"Unknown integration method {self.method!r}. Use one of: {valid}.")

        return self.sol

    def transform(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Convert angular coordinates to Cartesian positions."""
        if self.y is None:
            raise ValueError("Run the simulation first.")

        theta1, theta2 = self.y[0], self.y[1]

        x1 = self.params.L1 * np.sin(theta1)
        y1 = -self.params.L1 * np.cos(theta1)
        x2 = x1 + self.params.L2 * np.sin(theta2)
        y2 = y1 - self.params.L2 * np.cos(theta2)

        return x1, y1, x2, y2

    def energies(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return kinetic, potential, and total mechanical energy."""

        if self.y is None:
            raise ValueError("Run the simulation first.")

        theta1, theta2, omega1, omega2 = self.y
        m1, m2 = self.params.m1, self.params.m2
        L1, L2 = self.params.L1, self.params.L2
        g = self.params.g

        _, y1, _, y2 = self.transform()

        vx1 = L1 * omega1 * np.cos(theta1)
        vy1 = L1 * omega1 * np.sin(theta1)
        vx2 = vx1 + L2 * omega2 * np.cos(theta2)
        vy2 = vy1 + L2 * omega2 * np.sin(theta2)

        kinetic_1 = 0.5 * m1 * (vx1**2 + vy1**2)
        kinetic_2 = 0.5 * m2 * (vx2**2 + vy2**2)
        kinetic = kinetic_1 + kinetic_2

        potential_1 = m1 * g * y1
        potential_2 = m2 * g * y2
        potential = potential_1 + potential_2

        return kinetic, potential, kinetic + potential