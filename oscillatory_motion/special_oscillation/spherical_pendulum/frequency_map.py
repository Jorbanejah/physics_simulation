"""
Special graphics:

- Frequency-maps analysis or Frequency Laskar maps.

    - 1_ Write the equation 
    - 2_ Choose a grid of initial conditions: 
            1_ Choose an energy level
            2_ Scan over initial conditions angles
    - 3_Integrate the equation long enough times. (DOP853 with t_times = 250)
    - 4_Build a complex signal for frequency extraction -> E.g: Z (t) = theta(t) + j * phi(t)
    - 5_Apply Laskar-style frequency analysis (NAFF): use a Fourier to approximate Z(t), extract the fundamental omega_k and take a frequency vector (v_1, v_2)
    - 6_Construct the frequency map: associate the initial conditions with its frequencies, plot the frequencies vector v_1/ v_2vs the peak parameter with a fancy colormap


    A POWERFUL Laskar map could be: 
        
    - Two-window analysis: compute frequencies over two successive time windows for each trajectory.

    - Frequency drift: measure  Δv = v(1)-v(2); large drifts signal chaotic diffusion.


    Friendly note: the polar angle is phi and the axial angle is theta 
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, Any, Tuple, Callable
from collections.abc import Sequence

 

#Define the current state
State = np.ndarray

Dynamics = Callable[[float, State, "Params"], State] #Input -> [time, state, params]; Output -> [State]
##
# -------------------------- Params and quations -----------------------
##
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
    p0: tuple[float, float] = (0.0, 0.0)

    t:float = 250
    dt: float = 0.01

    def __post_init__(self) -> None:
        if self.g <= 0.0:
            raise ValueError("Gravity mus be positive.")
        if self.m <= 0.0:
            raise ValueError("Mass must be positive.")
        if self.L <= 0.0:
            raise ValueError("Lenghts must be positive.")
        
        self.q0 = _as_pair(self.q0, "q0")
        self.p0 = _as_pair(self.p0, "dq0")

def spherical_pendulum_equation(t: float, y:State, p:Params) -> State:
    """
    Full nonlinear differential spherical-pendulum equations (hamiltonian)

    State vector: [theta, phi, mtheta, mphi] where the mtheta and mphi correspond to the conjugate momentum
    define through hamiltonian equation.
    """

    del t

    m, L, g = p.m, p.L, p.g

    theta, phi, mtheta, mphi = np.asarray(y, dtype = float)

    eps = 1e-8
    if np.abs(phi) <eps:
        phi = np.sign(phi) * eps

    dtheta = mtheta/(m * L **2)
    dphi = mphi/(m *L**2 * np.sin(phi))

    dmphi = mphi**2/(m * L**2 * np.sin(phi) **3) *np.cos(phi) - m * g * L * np.sin(phi) 
    dmtheta = 0
   
    return np.array([dtheta, dphi, dmtheta, dmphi])


##
# ------------------- Energy level, initial conditions grid and integration --------
##

def energy(y: State, p: Params)-> float:

    theta, phi, mtheta, mphi = np.asarray(y, dtype = float)
    m, L, g = p.m, p.L, p.g
    T = mtheta**2/(2*m*L**2) + mphi**2/(2*m*L**2*np.sin(phi)**2)
    V = m*g*L*np.cos(phi)

    return T + V

# Knowing that the minimun energy is when phi = 0 -> E_min = mgL then my energy level have to be E_0 > E_min.
# After that, we are going to consider three kind of level energy:
# small level energy (E0 + deltaE, where deltaE << mgL), 
# medium level energy (E0 can reach between 30-60) 
# strong level energy (E0 can reach until pi/2)

def classify_energy(Et: float, params: Params) -> str:
    "The function classify the current energy of the state base on potential energy"
    m, L, g = params.m, params.L, params.g

    E_min = m*g*L

    # Energy thresholds based on maximum reachable polar angle.
    E_small  = m*g*L*np.cos(np.deg2rad(15))   # small oscillations
    E_medium = m*g*L*np.cos(np.deg2rad(45))   # medium oscillations
    E_strong = m*g*L*np.cos(np.deg2rad(90))   # strong oscillations

    if Et >= E_small:
        return "small"
    elif Et >= E_medium:
        return "medium"
    else:
        return "strong"

##
# ---------------------- Generate initial conditions ----------------
##
# The tricky part is, given a E0 energy, you have to find the initial conditions who generate this surface. 
# Given a theta_0, and phi_0, the unique form to reach the level E0 is through conjugate momenta

def generate_initial_conditions(E0: float, params: Params, N: int =50)-> np.array:
    """
    Description
    -------------
    The function calculate the family of initial conditions [theta, phi, mtheta, mphi] that describe the energy surface E0.
    How does it work? The program choose N theta and phi as initial angle, then, for simplicity, mtheta = 0 -conjugate momentum of theta.
    The following steps are:

        - Calculate the current potential energy: V = mgL cos(phi)
        - Calculate the T_needed = E0 - V
        - From kinetic energy equation, we solve for the conjugate momentum of phi as: 
        mphi = np.sqrt(2*m*L**2*np.sin(ph)**2 * T_needed)

    Parameters
    --------------
    E0: float
        Energy surface
    params: Params
        The current parameters which describe the motion
    N: int
        Grid configuration N x N
    """

    m, L, g = params.m, params.L, params.g

    thetas = np.linspace(0, 2*np.pi, N)
    phis   = np.linspace(0.05, np.pi-0.05, N)

    initials = []
  
    for th in thetas:
        for ph in phis:
           
            # Once you decide the initial angle, you have to choose the momenta to energy = E0
            # Choose mtheta = 0 for simplicity
            mtheta = 0.0

            # Solve for mphi from energy equation
            V = m*g*L*np.cos(ph)
            T_needed = E0 - V

            if T_needed <= 0:
                continue

            mphi = np.sqrt(2*m*L**2*np.sin(ph)**2 * T_needed)

            initials.append(np.array([th, ph, mtheta, mphi]))

    return initials

##
# ------------------- Integration -------------------
##

def integrate_trajectory(y0: State, params: Params) -> Tuple[Any, Any]:

    tgrid = np.arange(0, params.t, params.dt)
    sol = solve_ivp(lambda t, y: spherical_pendulum_equation(t, y, params), [0, params.t], y0, t_eval=tgrid, rtol=1e-9, atol=1e-9)

    return sol.t, sol.y

##
# ---------------- Complex signal and extract the frequencies -------------
##

def build_signal(theta: Sequence[float], phi: Sequence[float]):
    return theta + 1j*phi

def extract_frequency(signal, dt: float)-> Sequence[float]:
    """
    Through the Discrete Fourier transform we can extract the different frequencies from the current signal.
    This is not a complete NAFF - is too complex and I do not understand at all.
    """
    S = np.fft.fft(signal)
    freqs = np.fft.fftfreq(len(signal), dt)

    # Dominant frequency
    idx = np.argmax(np.abs(S))
    return abs(freqs[idx])

def compute_frequencies(t:Sequence[float], y:State) -> Tuple:

    "For a spherical pendulum, we need two angles so two frequencies. "

    theta, phi = y[0], y[1] #Extract the current trajectory

    #Build both complex signal
    sig_theta = build_signal(theta, phi)
    sig_phi   = build_signal(phi, theta)

    #Extract the current frequencies for each signal
    f1 = extract_frequency(sig_theta, t[1]-t[0])
    f2 = extract_frequency(sig_phi, t[1]-t[0])

    return f1, f2

##
# ------------------------- Construct the frequency map and plot it --------------
##

def frequency_map(initial_conditions: Sequence[float], params: Params):
    freqs = []
    count =0 
    total_steps = len(initial_conditions)
    for y0 in initial_conditions:
        count +=1
        print(f"Grid Progress: {100*count/total_steps:.1f}%", end="\r")
        t, y = integrate_trajectory(y0, params)
        f1, f2 = compute_frequencies(t, y)
        freqs.append((f1, f2))

    return np.array(freqs)

def plot_frequency_map(freqs):
    f1 = freqs[:,0]
    f2 = freqs[:,1]

    plt.figure(figsize=(8,6))
    plt.scatter(f1, f2, s=5, c='blue')
    plt.xlabel("Frequency 1")
    plt.ylabel("Frequency 2")
    plt.title("Laskar Frequency Map")
    plt.grid(True)
    plt.show()


##
# -------------------------- Main --------------------
##

params = Params()

# Choose energy level
E0 = params.m * params.g * params.L * np.cos(np.deg2rad(45))  # medium energy

# Generate initial conditions
print("Generating the initial conditions")
initials = generate_initial_conditions(E0, params, N=10)

# Compute frequency map
print("Computing the frequency map")
freqs = frequency_map(initials, params)

# Plot
print("Plotting")
plot_frequency_map(freqs)
