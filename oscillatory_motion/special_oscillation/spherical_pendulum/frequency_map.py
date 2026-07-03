"""

├── Params
├── spherical_pendulum_equations()
├── energy()
├── generate_initial_conditions()
├── integrate_trajectory()
│
├── hann_window()
├── analytic_signal()
├── fft_guess()
├── scalar_product()
├── refine_frequency()
├── naff()
├── naff_decomposition()
├── frequency_drift()
│
├── compute_frequency_map()
├── plot_frequency_map()
│
└── main()

    Friendly note: the polar angle is phi and the azimuthal angle is theta 
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.signal import hilbert
from scipy.optimize import minimize_scalar
from dataclasses import dataclass
from typing import Any, Tuple, Callable
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

    Friendly reminder: theta is the azimuthal angle and phi the polar angle
    """

    del t

    m, L, g = p.m, p.L, p.g

    theta, phi, mtheta, mphi = np.asarray(y, dtype = float)

    eps =1e-8 # Avoid the singularity

    sin_phi = np.sin(phi)

    if abs(sin_phi) < eps:
        sin_phi = eps if sin_phi >= 0 else -eps

    #Generalized velocities
    dphi = mphi/(m * L **2)
    dtheta = mtheta/(m *L**2 * sin_phi **2)

    dmphi = (mtheta**2 *np.cos(phi))/(m * L**2 * sin_phi **3) - m * g * L * np.sin(phi) 
    dmtheta = 0

    #if abs(dmphi) > 1e6:
    #    print(phi, sin_phi, dmphi)
   
    return np.array([dtheta, dphi, dmtheta, dmphi])


##
# ------------------- Energy, energy levels and initial conditions grid --------
##

def energy(y: State, p: Params)-> float:

    theta, phi, mtheta, mphi = np.asarray(y, dtype = float)

    m, L, g = p.m, p.L, p.g
    T = mphi**2/(2*m*L**2) + mtheta**2/(2*m*L**2*np.sin(phi)**2)
    V = - m*g*L*np.cos(phi)

    return T + V

# Knowing that the minimun energy is when phi = 0 -> E_min = -mgL then my energy level have to be E_0 > E_min.
# After that, we are going to consider three kind of level energy:
# small level energy (E0 + deltaE, where deltaE << mgL), 
# medium level energy (E0 can reach between 30-60) 
# strong level energy (E0 can reach until pi/2)

def classify_energy(Et: float, params: Params) -> str:
    "The function classify the current energy of the state base on potential energy"
    m, L, g = params.m, params.L, params.g

    E_min = m * g * L

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


# The tricky part is, given a E0 energy, you have to find the initial conditions who generate this surface. 
# Given a theta_0, and phi_0, the unique form to reach the level E0 is through conjugate momenta

def generate_initial_conditions(E0: float, params: Params, Ntheta: int =50, Nphi:int = 50)-> list[np.ndarray]:
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
    Ntheta and Nphi: int, int
        Grid configuration Ntheta x Nphi
    """

    m, L, g = params.m, params.L, params.g

    phi_max = np.arccos(-E0/(m*g*L)) # Threshold

    phis = np.linspace(0.05, phi_max, Nphi)
    thetas = np.linspace(-np.pi, np.pi, Ntheta)
    i =0
    initials = []
    for th in thetas:
        for ph in phis:
           
            # Once you decide the initial angle, you have to choose the momenta to energy = E0
            # Choose mtheta = 0 and th =0.0 for simplicity
            mtheta = 0.0

            # Solve for mphi from energy equation
            V = -m*g*L*np.cos(ph)
            T_needed = E0 - V

            if T_needed <= 1e-12:
                i +=1
                continue

            mphi1 = np.sqrt(2*m*L**2*np.sin(ph)**2 * T_needed)
            mphi2 = -np.sqrt(2*m*L**2*np.sin(ph)**2 * T_needed) #-> if you descomment this line you have to store in initials list the whole results too. In this way, the odds will be the -mphi and evens the mphi results

            initials.append(np.array([th, ph, mtheta, mphi1]))
            initials.append(np.array([th, ph, mtheta, mphi2]))
    print(f"Trajectories discard for non-positive kinetic energy: {i}")
    return initials

##
# ------------------- Trajectory integration -------------------
##
@dataclass
class Trajectory:
    t: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    mtheta: np.ndarray
    mphi: np.ndarray

# Angle wrapped
def unwrap_angle(theta, phi):
    theta = np.unwrap(theta)

    return theta, phi

def integrate_trajectory(y0: State, params: Params) -> Trajectory:

    t = _time_grid(duration= params.t, dt = params.dt)
  
    sol = solve_ivp(lambda t, y: spherical_pendulum_equation(t, y, params), [0, params.t], y0, t_eval=t, method = "DOP853", rtol=1e-9, atol=1e-9)
    
    if not sol.success:
        raise RuntimeError(sol.message)
    
    theta = sol.y[0]
    phi = sol.y[1]

    theta, phi = unwrap_angle(theta, phi)

    #Check the drift energy and rejects those trajectories whose drift energy exceeds a tolerance

    tol = 1e-2
    Et = energy(y = sol.y, p = params)

    if np.abs(max(Et) - min(Et)) > tol:
        raise ValueError("The energy drift exceeds the tolerance. This trajectory will be rejected")
    
    else:
        return Trajectory(sol.t, theta, phi, sol.y[2], sol.y[3])

def split_trajectory(traj: Trajectory):

    n = len(traj.t)

    mid = n//2

    first = Trajectory(traj.t[:mid], traj.theta[:mid], traj.phi[:mid], traj.mtheta[:mid], traj.mphi[:mid])

    second = Trajectory(traj.t[mid:], traj.theta[mid:], traj.phi[mid:], traj.mtheta[mid:], traj.mphi[mid:])
    
    return first, second
##
# ---------------- Complex canonical signal -------------
##

def canonical_signals(traj: Trajectory):

    """
    Canonical complex varaibles

    We describe two complex signal.
    """

    z_theta = traj.theta  -1j *traj.mtheta
    z_phi = traj.phi -  1j * traj.mphi

    return z_theta, z_phi

##
# ---------------------- NAFF algorithms ----------------------
# - Window function: we use a Hann window. However, Laskar typically uses a high-order cosine window
# (order 4 or even 6). This improves convergence from approximately O(T-2)
#
# - After extracting each frequency, the original NAFF algorithm re-orthogonalizes the basis and recomputes all previosly found amplitudes.
# Simply subtracting Ae^(iwt) is good approximation but not mathematically exact.
#
# - Cosine window
# - fft_guess
# - scalar product
# - Refine frequencies
#
#
##

def cosine_window(N: int, order: int = 4):

    """
    Generalized cosine window used in NAFF.

    Parameters
    ----------
    order : int
        Window order.
        order=1 -> Hann window
        order=4 -> typical Laskar choice
    """

    x = np.linspace(-1.0, 1.0, N)

    window = (1.0 + np.cos(np.pi*x))**order

    # normalize so <1,1>=1

    window /= np.sum(window)

    return window

def fft_guess(signal, t):
    "FFT is only used to produce an initial guess"

    dt = t[1] - t[0]

    signal = signal - np.mean(signal)

    signal *= cosine_window(len(signal))

    F = np.fft.fft(signal)

    freq = np.fft.fftfreq(len(signal), dt)

    positive = freq > 0

    spectrum = F[positive]

    freq = freq[positive]

    k = np.argmax(np.abs(spectrum))

    return 2*np.pi*freq[k]

def scalar_product(signal, omega, t, window):
    """
    Computes the dot product using the Hann-weighted scalar product
    """
    exponential = np.exp(-1j*omega*t)

    return np.vdot(window*exponential, signal)

def refine_frequency(signal, omega0, t):

    w = cosine_window(len(signal))

    domega = 2*np.pi/(t[-1]-t[0])

    result = minimize_scalar(
        lambda om: np.abs(scalar_product(signal, om, t, w)),
        bounds=(omega0-domega, omega0+domega),
        method="Bounded"
    )

    return result.x

# ----------------- Refine amplitude ------------------

def estimate_amplitude(signal, omega, t, window, ):
    """
    Projection of the signal over e^(iwt).

    Returns the complex Fourier coefficient.
    """
    exponential = np.exp(-1j*omega*t)

    numerator = np.vdot(exponential*window, signal)

    denominator = np.vdot(exponential*window, exponential)
    
    return numerator/denominator

# ---------------------- One-frequency NAFF -------------------

def naff(signal, t,):
    """
    Extract the dominant frequency

    Returns (omega, amplitude)
    """
    signal = signal.astype(complex)

    signal -= np.mean(signal)

    window = cosine_window(len(signal), order =4)

    omega0 = fft_guess(signal, t)

    omega = refine_frequency(signal, omega0, t)

    amplitude = estimate_amplitude(signal, omega, t, window)

    return omega, amplitude

# -------------------- Multi-frequency extraction -------------

def naff_decomposition(signal, t, nfreq=5):
    """
    Iteratively extracts the dominant frequencies.
    """
    residual = signal.copy().astype(complex)

    omegas = []

    amplitudes = []

    for _ in range(nfreq):

        omega, A = naff(residual, t)

        residual -= A*np.exp(1j*omega*t)

        omegas.append(omega)

        amplitudes.append(A)

    return np.asarray(omegas), np.asarray(amplitudes)

# ------------------ Fundamental frequencies and frequencies diffusion ----------------

def fundamental_frequencies(trajectory: Trajectory,):

    """
    Computes the two fundamental frequencies of the spherical pendulum.

    Returns: omega_theta, omega_phi
    """

    z_theta, z_phi = canonical_signals(trajectory)

    omega_theta, _ = naff(z_theta, trajectory.t)

    omega_phi, _ = naff(z_phi, trajectory.t)

    return (omega_theta, omega_phi)

def frequency_drift(trajectory: Trajectory, ):

    """
    Computes Laskar's chaos indicator. 
    Delta = abs(omega1 - omega2) / omega1
    """

    first, second = split_trajectory(trajectory)

    w1_theta, w1_phi = fundamental_frequencies(first)

    w2_theta, w2_phi = fundamental_frequencies(second)

    drift_theta = np.abs(w2_theta - w1_theta)/np.abs(w1_theta)
    drift_phi = np.abs(w2_phi - w1_phi)/np.abs(w1_phi)

    return (w1_theta, w1_phi, drift_theta, drift_phi)
##
# ------------------------- Construct the frequency map and plot it --------------
##

@dataclass(slots=True)

class FrequencyPoint:

    theta0: float
    phi0: float

    omega_theta: float
    omega_phi: float

    drift_theta: float
    drift_phi: float

# Analyse one trajectory

def analyse_trajectory(y0, params,):
    """
    Performs the complete frequency analysis
    of one trajectory.
    """

    traj = integrate_trajectory(y0,params)

    omega_theta, omega_phi, drift_theta, drift_phi = frequency_drift(traj)

    return FrequencyPoint(theta0=y0[0], phi0=y0[1],
        omega_theta=omega_theta, omega_phi=omega_phi, 
        drift_theta=drift_theta, drift_phi=drift_phi,)

# Entire grid

def compute_frequency_map(initial_conditions, params):
    """
    Computes the frequency map for every initial condition.
    """
    results = []

    total = len(initial_conditions)

    for k, y0 in enumerate(initial_conditions):

        print(f"\r{k+1}/{total}", end="", flush=True,)

        try:

            point = analyse_trajectory(y0,params)

            results.append(point)

        except Exception as err:

            print()

            print("Skipped trajectory:", err)

    print()

    return results

# Convert to arrays

def unpack_results(results):
    theta0 = np.array(
        [r.theta0 for r in results])
    
    phi0 = np.array(
        [r.phi0 for r in results])
    
    omega_theta = np.array(
        [r.omega_theta for r in results])
    
    omega_phi = np.array(
        [r.omega_phi for r in results])
    
    drift_theta = np.array(
        [r.drift_theta for r in results])
    
    drift_phi = np.array(
        [r.drift_phi for r in results])
    
    return (theta0, phi0, omega_theta, omega_phi, drift_theta, drift_phi,)

# Frequency-frequency map

def plot_frequency_map(results):

    theta0, phi0, omega_theta, omega_phi, drift_theta, drift_phi, = unpack_results(results)
    ratio = omega_theta / omega_phi
    plt.figure(figsize=(8,7))
    plt.scatter(omega_theta, omega_phi, c=ratio, cmap="turbo", s=18, label = r"$Ratio = \omega_\theta/\omega_\phi$")
    plt.xlabel(r"$\omega_\theta$")
    plt.ylabel(r"$\omega_\phi$")
    plt.title("Frequency Map")
    cbar = plt.colorbar()
    cbar.set_label(r"$\omega_\theta/\omega_\phi$")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Frequency diffusion map

def plot_diffusion_map(results):

    theta0, phi0, omega_theta, omega_phi, drift_theta, drift_phi, = unpack_results(results)
    diffusion = np.maximum(drift_theta, drift_phi,)
    diffusion = np.log10(diffusion + 1e-16)
    plt.figure(figsize=(8,7))
    plt.scatter(theta0, phi0, c=diffusion, cmap="inferno", s=20,)
    plt.xlabel(r"$\theta_0$")
    plt.ylabel(r"$\phi_0$")
    plt.title("Laskar Diffusion Map")
    cbar = plt.colorbar()
    cbar.set_label(r"$\log_{10}(\Delta\omega)$")
    plt.tight_layout()
    plt.show()

# Resonance map

def plot_resonance_map(results):

    theta0, phi0, omega_theta, omega_phi, drift_theta, drift_phi, = unpack_results(results)
    ratio = omega_theta / omega_phi
    plt.figure(figsize=(8,7))
    plt.scatter(theta0, phi0, c=ratio, cmap= "twilight", s =20)
    plt.xlabel(r"$\theta_0$")
    plt.ylabel(r"$\phi_0$")
    plt.title("Frequency Ratio")
    cbar = plt.colorbar()
    cbar.set_label(r"$\omega_\theta/\omega_\phi$")
    plt.tight_layout()
    plt.show()

# Frequency drift histogram

def plot_drift_histogram(results):
    theta0, phi0, omega_theta, omega_phi, drift_theta, drift_phi, = unpack_results(results)
    diffusion = np.maximum(drift_theta, drift_phi,)
    plt.figure(figsize=(7,5))
    plt.hist(np.log10(diffusion+1e-16),bins=40,)
    plt.xlabel(r"$\log_{10}(\Delta\omega)$")
    plt.ylabel("Counts")
    plt.title("Frequency Drift Distribution")
    plt.tight_layout()
    plt.show()

##
# -------------------------- Save results and load CSV --------------------
##

from pathlib import Path
import pandas as pd

# Save results

def results_to_dataframe(results):

    rows = []

    for r in results:

        rows.append({
            "theta0": r.theta0,
            "phi0": r.phi0,
            "omega_theta": r.omega_theta,
            "omega_phi": r.omega_phi,
            "drift_theta": r.drift_theta,
            "drift_phi": r.drift_phi,
        })

    return pd.DataFrame(rows)

# Save CSV

def save_results(results, filename):

    df = results_to_dataframe(results)

    df.to_csv(filename, index=False)

    print(f"Saved {len(df)} trajectories")

# Load CSV

def load_results(filename):

    return pd.read_csv(filename)

# Plot directly from dataframe

def plot_dataframe(df):

    diffusion = np.maximum(
        df["drift_theta"],
        df["drift_phi"]
    )

    diffusion = np.log10(diffusion + 1e-16)

    plt.figure(figsize=(8,7))
    plt.scatter(df["theta0"], df["phi0"], c=diffusion, cmap="inferno", s=15,)
    plt.xlabel(r"$\theta_0$")
    plt.ylabel(r"$\phi_0$")
    plt.title("Laskar Diffusion Map")
    cbar = plt.colorbar()
    cbar.set_label(r"$\log_{10}\Delta\omega$")
    plt.tight_layout()
    plt.show()

# 
# ------------------------ Main ----------------
# 

def main():

    params = Params(g=9.81, m=1.0, L=2.0, t=250, dt=0.01)

    # Choose the energy surface
    energy_level = (-params.m* params.g* params.L * np.cos(np.deg2rad(30)))

    # Initial-condition grid
    print("=" *60)

    initials = generate_initial_conditions(E0=energy_level, params=params, Nphi= 30, Ntheta=30)
    
    print()

    print("-"*20)
    print("Energy level")
    print(energy_level)
    print("Trajectories")
    print(len(initials))
    print("-"*20)

    # Frequency map
    results = compute_frequency_map(initials,params,)

    # Save    
    output = Path("frequency_map.csv")

    save_results(results, output,)


    plot_frequency_map(results)
    plot_diffusion_map(results)

    plot_resonance_map(results)
    plot_drift_histogram(results)

if __name__ == "__main__":

    main()