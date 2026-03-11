''''
1️⃣ trayectorias x(t) para los 3 regímenes
2️⃣ espacio de fase
3️⃣ energía mecánica
4️⃣ energía total
5️⃣ envolvente exponencial

Numéricas

6️⃣ error vs tiempo
7️⃣ error vs dt
8️⃣ error energético
9️⃣ error vs β
🔟 tiempo de relajación vs parámetro
'''

import numpy as np
import matplotlib.pyplot as plt
from Damping_vibration import Pendulum
from dataclasses import dataclass

@dataclass
class DampingPendulum:
    'Physical and numerical parameters'

    mass: float = 2
    theta0: float = np.deg2rad(15)
    omega0: float = 0
    gamma: float = 2
    dt = 0.01
    t_max: float = 15

# Only will be use for approximation
def analytic_solution(t, y0, v0, beta, omega0):
    """
    Analytical solution of the damped harmonic oscillator.
    """

    if beta < omega0:

        omega_d = np.sqrt(omega0**2 - beta**2)

        A = y0
        B = (v0 + beta*y0) / omega_d

        return np.exp(-beta*t) * (A*np.cos(omega_d*t) + B*np.sin(omega_d*t))

    elif beta == omega0:

        return np.exp(-beta*t) * (y0 + (v0 + beta*y0)*t)

    else:

        r1 = -beta + np.sqrt(beta**2 - omega0**2)
        r2 = -beta - np.sqrt(beta**2 - omega0**2)

        C1 = (v0 - r2*y0)/(r1 - r2)
        C2 = y0 - C1

        return C1*np.exp(r1*t) + C2*np.exp(r2*t)


def simulate_pendulum(pendulum, params, method):
    """
    Simulate the damped pendulum using a chosen numerical method.

    """

    theta = pendulum.theta
    omega = pendulum.omega

    times = [0]
    theta_vals = [theta]
    omega_vals = [omega]

    Ek, Ep, Wp = [], [], []

    Ek0, Ep0, W0, _ = pendulum.energy(theta, omega, 0)

    Ek.append(Ek0)
    Ep.append(Ep0)
    Wp.append(W0)

    relaxation_time = None

    for t in np.arange(params.dt, params.t_max, params.dt):

        # select integration method
        if method == "rk4":
            theta, omega = pendulum.rk4(params.dt)

        elif method == "euler":
            theta, omega = pendulum.euler(params.dt)

        elif method == "crank_nicolson":
            theta, omega = pendulum.crank_nicolson(params.dt)

        else:
            raise ValueError("Unknown method")

        times.append(t)
        theta_vals.append(theta)
        omega_vals.append(omega)

        Ek0, Ep0, W0, _ = pendulum.energy(theta, omega, t)

        Ek.append(Ek0)
        Ep.append(Ep0)
        Wp.append(Wp[-1] + W0)

        # relaxation condition
        if relaxation_time is None:

            if abs(theta) <= abs(params.theta0) * np.exp(-1):

                relaxation_time = t

    Et = np.array(Ek) + np.array(Ep) + np.array(Wp)

    return {
        "times": np.array(times),
        "theta": np.array(theta_vals),
        "omega": np.array(omega_vals),
        "Ek": np.array(Ek),
        "Ep": np.array(Ep),
        "Wp": np.array(Wp),
        "Et": np.array(Et),
        "relaxation": relaxation_time
    }

##
# Plots
##


# ---------------------------------------------------------
# 1. TRAJECTORY
# ---------------------------------------------------------
def plot_trajectory(times, x, label=None):

    plt.figure(figsize=(8,5))

    plt.plot(times, x, lw=2, label=label)

    plt.xlabel("Time (s)")
    plt.ylabel("x(t)")
    plt.title("Trajectory")
    plt.grid(alpha=0.3)

    if label:
        plt.legend()

    plt.tight_layout()


# ---------------------------------------------------------
# 2. PHASE SPACE
# ---------------------------------------------------------
def plot_phase_space(x, v, label=None):

    plt.figure(figsize=(6,6))

    plt.plot(x, v, lw=2, label=label)

    plt.xlabel("x")
    plt.ylabel("v")
    plt.title("Phase space")
    plt.grid(alpha=0.3)

    if label:
        plt.legend()

    plt.tight_layout()



# ---------------------------------------------------------
# 3. ENERGY
# ---------------------------------------------------------
def plot_energy(times, Ek, Ep, Et=None):

    plt.figure(figsize=(8,5))

    plt.plot(times, Ek, label="Kinetic energy")
    plt.plot(times, Ep, label="Potential energy")

    if Et is not None:
        plt.plot(times, Et, lw=2, label="Total energy")

    plt.xlabel("Time (s)")
    plt.ylabel("Energy")
    plt.title("Energy evolution")

    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()


# ---------------------------------------------------------
# 4. NUMERICAL ERROR vs TIME
# ---------------------------------------------------------
def plot_error_time(times, x_num, x_exact):

    error = np.abs(x_num - x_exact)

    plt.figure(figsize=(8,5))

    plt.plot(times, error)

    plt.xlabel("Time (s)")
    plt.ylabel("Error |x_num - x_exact|")
    plt.title("Numerical error")

    plt.yscale("log")
    plt.grid(alpha=0.3)
    plt.tight_layout()


# ---------------------------------------------------------
# 5. ERROR vs dt
# ---------------------------------------------------------
def plot_error_dt(dt_values, errors):

    plt.figure(figsize=(7,5))

    plt.loglog(dt_values, errors, marker="o")

    plt.xlabel("dt")
    plt.ylabel("Error")
    plt.title("Convergence test")

    plt.grid(True, which="both")
    plt.tight_layout()


# ---------------------------------------------------------
# 6. RELAXATION TIME
# ---------------------------------------------------------
def plot_relaxation(beta, tau):

    plt.figure(figsize=(7,5))

    plt.plot(beta, tau, marker="o")

    plt.xlabel("Damping coefficient β")
    plt.ylabel("Relaxation time τ")

    plt.title("Relaxation time vs damping")

    plt.grid(alpha=0.3)
    plt.tight_layout()


# ---------------------------------------------------------
# 7. ENERGY ERROR
# ---------------------------------------------------------
def plot_energy_error(times, energy):

    error = np.abs(energy - energy[0])

    plt.figure(figsize=(8,5))

    plt.plot(times, error)

    plt.xlabel("Time (s)")
    plt.ylabel("|E(t) - E(0)|")

    plt.title("Energy drift")

    plt.yscale("log")
    plt.grid(alpha=0.3)
    plt.tight_layout()


# ---------------------------------------------------------
# 8. THREE REGIMES COMPARISON
# ---------------------------------------------------------
def plot_three_regimes(times_list, x_list, labels):

    plt.figure(figsize=(9,6))

    for t, x, lab in zip(times_list, x_list, labels):

        plt.plot(t, x, lw=2, label=lab)

    plt.xlabel("Time (s)")
    plt.ylabel("x(t)")
    plt.title("Damping regimes comparison")

    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()


# ---------------------------------------------------------
# 9. PHASE SPACE THREE REGIMES
# ---------------------------------------------------------
def plot_phase_three(x_list, v_list, labels):

    plt.figure(figsize=(6,6))

    for x, v, lab in zip(x_list, v_list, labels):

        plt.plot(x, v, label=lab)

    plt.xlabel("x")
    plt.ylabel("v")
    plt.title("Phase space comparison")

    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()