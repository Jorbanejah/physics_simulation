import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from Damping_vibration import Pendulum


# ---------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------

@dataclass
class DampingPendulum:
    mass: float = 2
    gamma: float = 1
    L: float = 1
    g: float = 9.81

    theta0: float = np.deg2rad(10)
    omega0: float = 0

    dt: float = 0.01
    t_max: float = 15


# ---------------------------------------------------------
# ANALYTIC SOLUTION (SMALL ANGLE)
# ---------------------------------------------------------

def analytic_solution(t, theta0, omega0, beta, w0):

    if beta < w0:

        wd = np.sqrt(w0**2 - beta**2)

        A = theta0
        B = (omega0 + beta*theta0)/wd

        return np.exp(-beta*t)*(A*np.cos(wd*t) + B*np.sin(wd*t))

    elif beta == w0:

        return np.exp(-beta*t)*(theta0 + (omega0+beta*theta0)*t)

    else:

        r1 = -beta + np.sqrt(beta**2 - w0**2)
        r2 = -beta - np.sqrt(beta**2 - w0**2)

        C1 = (omega0-r2*theta0)/(r1-r2)
        C2 = theta0 - C1

        return C1*np.exp(r1*t) + C2*np.exp(r2*t)


# ---------------------------------------------------------
# SIMULATION
# ---------------------------------------------------------

def simulate_pendulum(pendulum, params, method):

    theta = pendulum.theta
    omega = pendulum.omega

    times=[0]
    theta_vals=[theta]
    omega_vals=[omega]

    Ek=[]
    Ep=[]
    Wp=[]

    Ek0,Ep0,W0,_=pendulum.energy(theta,omega,0)

    Ek.append(Ek0)
    Ep.append(Ep0)
    Wp.append(W0)

    relaxation=None

    for t in np.arange(params.dt, params.t_max, params.dt):

        if method=="rk4":

            theta, omega=pendulum.rk4(params.dt)

        elif method=="euler":

            theta ,omega=pendulum.euler(params.dt)

        elif method=="crank_nicolson":

            theta, omega=pendulum.crank_nicolson(params.dt)

        else:

            raise ValueError("Unknown method")

        times.append(t)
        theta_vals.append(theta)
        omega_vals.append(omega)

        Ek0,Ep0,W0,_=pendulum.energy(theta,omega,params.dt)

        Ek.append(Ek0)
        Ep.append(Ep0)
        Wp.append(Wp[-1] + W0)

        pendulum.theta = theta
        pendulum.omega = omega

        if relaxation is None:

            if abs(theta)<=abs(params.theta0)*np.exp(-1):

                relaxation=t

    Et=np.array(Ek) + np.array(Ep) + np.array(Wp)

    return dict(
        times=np.array(times),
        theta=np.array(theta_vals),
        omega=np.array(omega_vals),
        Ek=np.array(Ek),
        Ep=np.array(Ep),
        Wp=np.array(Wp),
        Et=Et,
        relaxation=relaxation
    )


# ---------------------------------------------------------
# MAIN ANALYSIS
# ---------------------------------------------------------

def main_physics():

    params=DampingPendulum()

    methods=["rk4","euler","crank_nicolson"]

    gamma_values=[0.2 ,2 ,6]

    results={}

    for gamma in gamma_values:

        params.gamma=gamma

        results[gamma]={}

        for method in methods:

            pend=Pendulum(
                params.theta0,
                params.omega0,
                params.mass,
                params.L,
                params.gamma,
                params.dt,
                params.t_max,
                approx=False
            )

            results[gamma][method]=simulate_pendulum(
                pend,
                params,
                method
            )
    plot_regime_summary(results, gamma_values, methods)
    plt.show()

def main_errors():

    """
    Numerical error analysis for the damped pendulum.

    Includes:

    1. Energy conservation vs initial angle
    2. Convergence test using RK4 reference solution
    3. Numerical stability test
    """

    params = DampingPendulum()

    methods = ["rk4", "euler", "crank_nicolson"]

    dt_values = [0.1, 0.05, 0.01]

    gamma_values = [0.2, 2.0, 6.0]

    theta_deg_values = np.arange(1, 51)
    theta_values = np.deg2rad(theta_deg_values)

    params.dt = 0.01

    results_error = {}

# -----------------------------------------------------
# SIMULATION LOOP
# -----------------------------------------------------

    for gamma in gamma_values:

        params.gamma = gamma
        results_error[gamma] = {}

        for theta_deg, theta0 in zip(theta_deg_values, theta_values):

            results_error[gamma][theta_deg] = {}

            for method in methods:

                params.theta0 = theta0

                pend = Pendulum(
                    params.theta0,
                    params.omega0,
                    params.mass,
                    params.L,
                    params.gamma,
                    params.dt,
                    params.t_max,
                    approx=False
                )

                results_error[gamma][theta_deg][method] = simulate_pendulum(
                    pend,
                    params,
                    method
                )

# -----------------------------------------------------
# ENERGY CONSERVATION ANALYSIS
# -----------------------------------------------------

    plot_energy_vs_initial_angle(
        theta_deg_values,
        results_error,
        methods,
        gamma_values
    )

# -----------------------------------------------------
# CONVERGENCE ANALYSIS (REFERENCE RK4)
# -----------------------------------------------------

    theta_test = np.deg2rad(10)

    dt_ref = 1e-4

    params.dt = dt_ref
    params.theta0 = theta_test

    pend_ref = Pendulum(
        params.theta0,
        params.omega0,
        params.mass,
        params.L,
        params.gamma,
        params.dt,
        params.t_max,
        approx=False
    )

    reference = simulate_pendulum(pend_ref, params, "rk4")

    theta_ref = reference["theta"]
    t_ref = reference["times"]

    for method in ["rk4", "crank_nicolson"]:

        errors = []

        for dt in dt_values:

            params.dt = dt
            params.theta0 = theta_test

            pend = Pendulum(
                params.theta0,
                params.omega0,
                params.mass,
                params.L,
                params.gamma,
                params.dt,
                params.t_max,
                approx=False
            )

            r = simulate_pendulum(pend, params, method)

            theta_num = r["theta"]
            t_num = r["times"]

            # interpolate reference solution
            theta_ref_interp = np.interp(t_num, t_ref, theta_ref)

            error = np.sqrt(np.mean((theta_num - theta_ref_interp) ** 2))

            errors.append(error)

        plot_convergence(dt_values, errors, method)

# -----------------------------------------------------
# NUMERICAL STABILITY
# -----------------------------------------------------

    theta_test = np.deg2rad(10)

    for method in ["rk4", "crank_nicolson"]:

        amplitudes = []

        for dt in dt_values:

            params.dt = dt
            params.theta0 = theta_test

            pend = Pendulum(
                params.theta0,
                params.omega0,
                params.mass,
                params.L,
                params.gamma,
                params.dt,
                params.t_max,
                approx=False
            )

            r = simulate_pendulum(pend, params, method)

            amplitudes.append(np.max(np.abs(r["theta"])))

        plot_stability(dt_values, amplitudes, method)

    plt.show()

def plot_regime_summary(results, gamma_values, methods):

    """
    Create a 3x3 subplot showing the dynamics of the damped pendulum
    for the three damping regimes.

    Rows:
        1 -> trajectory θ(t)
        2 -> phase space (θ, ω)
        3 -> total energy E(t)

    Columns:
        underdamped | critical | overdamped
    """

    fig, axes = plt.subplots(3, 3, figsize=(15, 10))

    regime_names = [
        "Underdamped",
        "Critical damping",
        "Overdamped"
    ]
    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(methods)* 2) for i in range(len(methods))] # Parameter five is used to create a certain distance between the chosen colours.

    for col, gamma in enumerate(gamma_values):

        for i, method in enumerate(methods):

            r = results[gamma][method]
            # -------------------------
            # Row 1 : Trajectory
            # -------------------------
            axes[0, col].plot(r["times"], r["theta"], color = colors[i], label=method)

            # -------------------------
            # Row 2 : Phase space
            # -------------------------
            axes[1, col].plot( r["theta"],r["omega"], color = colors[i], label=method)

            # -------------------------
            # Row 3 : Energy
            # -------------------------
            axes[2, col].plot( r["times"],r["Et"], color = colors[i], label=method)

        # Column titles
        axes[0, col].set_title(regime_names[col] + rf'$\gamma ={gamma}$')

        # axis labels
        axes[0, col].set_xlabel("Time (s)")
        axes[0, col].set_ylabel("Angle (rad)")

        axes[1, col].set_xlabel("Angle (rad)")
        axes[1, col].set_ylabel("Angular velocity (rad/s)")

        axes[2, col].set_xlabel("Time (s)")
        axes[2, col].set_ylabel("Total energy (J)")

        for row in range(3):
            axes[row, col].grid(alpha=0.3)

    # legend only once
    axes[0,2].legend()

    plt.tight_layout()

def plot_energy_vs_initial_angle(theta_deg_values, results_error, methods, gamma_values):
    """
    Plot maximum and minimum total energy vs initial angle.

    A 1x3 subplot is created, one for each damping regime
    (different gamma values).
    """

    fig, ax = plt.subplots(3, 3, figsize=(15,5))

    regime_names = [
        "Underdamped",
        "Critical damping",
        "Overdamped"
    ]

    i = 0

    for gamma, regime in zip(gamma_values, regime_names):
        for j, method in enumerate(methods):

            max_energy = []
            min_energy = []

            for theta_deg in theta_deg_values:

                r = results_error[gamma][theta_deg][method]

                max_energy.append(np.max(r["Et"]))
                min_energy.append(np.min(r["Et"]))

            ax[i, j].plot(theta_deg_values, max_energy,
                    label=f"{method} max(E)")

            ax[i, j].plot(theta_deg_values, min_energy,
                    linestyle="--",
                    label=f"{method} min(E)")

            ax[i, j].set_title(rf"{regime} $( \gamma = {gamma})$")
            ax[i, j].set_xlabel(rf"Initial angle $\theta_0$  (deg)")
            ax[i ,j].grid(alpha=0.3)

            ax[i,j].set_ylabel("Total energy (J)")
            ax[i ,j].legend()

        i += 1

       

    fig.suptitle("Energy conservation vs initial angle")

    fig.tight_layout()


def plot_convergence(dt_values, errors, method):
    """
    Plot log(error) vs log(dt) to estimate the order of convergence.
    """

    plt.figure(figsize=(6,5))

    plt.loglog(dt_values, errors, "o-", label=method)

    plt.xlabel("Time step dt")
    plt.ylabel("Error")

    plt.title(f"Convergence test ({method})")

    plt.grid(True, which="both")
    plt.legend()
    plt.tight_layout()

def plot_stability(dt_values, amplitudes, method):
    """
    Plot numerical stability by measuring the maximum
    oscillation amplitude for different time steps.
    """

    plt.figure(figsize=(6,5))

    plt.plot(dt_values, amplitudes, "o-", label=method)

    plt.xlabel("Time step dt")
    plt.ylabel("Max |θ|")

    plt.title(f"Numerical stability ({method})")

    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()

if __name__=="__main__":
    main_physics()
    main_errors()