import numpy as np
import matplotlib.pyplot as plt
from time import time
from Driven_oscillation import DrivenOscillation


# ============================================================
#   Parameter container
# ============================================================

class DrivenOscillationParams:
    """
    Parameter set for a driven (linear or nonlinear) oscillator.
    """
    def __init__(self,
                 mass=1.0,
                 gamma=1.0,
                 L=1.0,
                 omega=2.0,
                 F0=2.0,
                 F_external='cos',
                 q0=np.deg2rad(30),
                 dq0=0.0,
                 dt=0.01,
                 t_max=50.0,
                 system='nonlinear'):

        self.mass = mass
        self.gamma = gamma
        self.L = L
        self.omega = omega
        self.F0 = F0
        self.F_external = F_external
        self.q0 = q0
        self.dq0 = dq0
        self.dt = dt
        self.t_max = t_max
        self.system = system

    def as_key(self, method: str):
        """
        Build a hashable key for caching.
        """
        return (self.q0, self.dq0,
                self.mass, self.gamma, self.L,
                self.F0, self.omega,
                self.dt, self.t_max,
                self.system, self.F_external,
                method)


# ============================================================
#   Simulation + cache
# ============================================================

_SIM_CACHE = {}


def run_simulation(params: DrivenOscillationParams):
    """
    Run a driven oscillation simulation and return (history, runtime).
    """
    start = time()

    osc = DrivenOscillation(q0 = params.q0, dq0 = params.dq0, m =params.mass, gamma=params.gamma, F0 = params.F0, omega= params.omega, t= params.t_max, dt = params.dt, system =params.system,L = params.L,F_external =params.F_external)

    model = osc.run()
    runtime = time() - start
    return model.history, runtime


def cached_simulation(params: DrivenOscillationParams, method: str):
    """
    Cached version of run_simulation to avoid recomputing identical cases.
    """
    key = params.as_key(method)
    if key in _SIM_CACHE:
        return _SIM_CACHE[key]

    history, runtime = run_simulation(params)
    _SIM_CACHE[key] = (history, runtime)
    return history, runtime


# ============================================================
#   Energy + diagnostics from history
# ============================================================

def total_energy(history, method: str):
    """
    Compute total energy: mechanical + work.
    """
    Ek = np.array(history[method]["Ek"])
    Ep = np.array(history[method]["Ep"])
    Wd = np.array(history[method]["Wp_diss"])
    Wf = np.array(history[method]["Wp_drive"])
    return Ek + Ep + (Wd - Wf)


def energy_drift_from_history(history, method: str):
    """
    Energy drift = max(E) - min(E) from a given history.
    """
    E = total_energy(history, method)
    return np.max(E) - np.min(E)


def max_amplitude_from_history(history, method: str):
    """
    Maximum absolute angle from a given history.
    """
    q = history[method]["q"]
    return np.max(np.abs(q))


# ============================================================
#   Plotting utilities
# ============================================================

def plot_energy_vs_initial_angle(theta, drift):
    plt.figure(figsize=(6, 5))
    plt.plot(theta, drift)
    plt.xlabel("Initial angle (rad)")
    plt.ylabel(r"Energy drift $\Delta E$")
    plt.title("Energy Drift vs Initial Angle")
    plt.grid(alpha=0.3)
    plt.tight_layout()

def plot_runtime(dt_values, runtimes):
    plt.figure(figsize=(6, 5))
    plt.plot(dt_values, runtimes, "o-", lw=2)
    plt.xlabel(r"Time step $\Delta t$")
    plt.ylabel("Runtime (s)")
    plt.title("Runtime vs Time Step")
    plt.grid(alpha=0.3)
    plt.tight_layout()

def plot_stability(dt_values, amplitudes, method):
    plt.figure(figsize=(6, 5))
    plt.plot(dt_values, amplitudes, "o-", lw=2, markersize=6)
    plt.xlabel(r"Time step $\Delta t$")
    plt.ylabel(r"Max $|\theta|$")
    plt.title(f"Numerical Stability ({method})", loc="left")
    plt.grid(alpha=0.3)
    plt.tight_layout()


def plot_convergence(dt_values, errors, method):
    """
    Plot log(error) vs log(dt) and include theoretical slope lines.
    """
    dt_values = np.array(dt_values, dtype=float)
    errors = np.array(errors, dtype=float)

    plt.figure(figsize=(6, 5))
    plt.loglog(dt_values, errors, "o-", lw=2, markersize=6, label=method)

    theoretical_order = {
        "rk4": 4,
        "crank_nicolson": 2,
        "verlet": 2
    }.get(method.lower(), None)

    if theoretical_order is not None:
        x0 = dt_values[len(dt_values) // 2]
        y0 = errors[len(errors) // 2]
        slope_line = y0 * (dt_values / x0) ** theoretical_order

        plt.loglog(dt_values, slope_line, "--", color="gray",
                   label=rf"Slope ≈ {theoretical_order}")

    plt.xlabel(r"Time step $\Delta t$")
    plt.ylabel("Error")
    plt.title(f"Convergence Test ({method})", loc="left")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()


# ============================================================
#   Main error analysis
# ============================================================

def run_error_analysis(base_method):
    """
    Run all numerical analyses for the driven nonlinear oscillator.

    Produces:
        1. Energy drift vs initial angle
        2. Convergence test
        3. Numerical stability test
        4. Runtime vs dt
    """

    params = DrivenOscillationParams()
    dt_values = [0.2, 0.1, 0.05, 0.01]

    results = {
        "energy_drift": {},
        "convergence": {},
        "stability": {},
        "runtime": {},
        "parameter_sweep": {}
    }

    # ---------------------------------------------------------
    # 1. Energy drift vs initial angle
    # ---------------------------------------------------------
    print("Energy drift vs initial angle")
    theta_deg = np.arange(1, 81)
    theta_values = np.deg2rad(theta_deg)

    drift_list = []
    for th in theta_values:
        params.q0 = th
        history, _ = cached_simulation(params, base_method)
        drift = energy_drift_from_history(history, base_method)
        drift_list.append(drift)

    drift_list = np.array(drift_list)
    results["energy_drift"]["theta"] = theta_values
    results["energy_drift"]["drift"] = drift_list
    plot_energy_vs_initial_angle(theta_values, drift_list)
    plt.show()

    # ---------------------------------------------------------
    # 2. Convergence test (reference RK4)
    # ---------------------------------------------------------
    print("Method convergence")
    errors = []
    for dt in dt_values:
        params.dt = dt
        history, _ = cached_simulation(params, base_method)
        drift = energy_drift_from_history(history, base_method)
        errors.append(drift)

    errors = np.array(errors)
    results["convergence"]["dt"] = np.array(dt_values)
    results["convergence"]["errors"] = errors
    plot_convergence(dt_values, errors, base_method)
    plt.show()
    # ---------------------------------------------------------
    # 3. Numerical stability test (max |theta| vs dt)
    # ---------------------------------------------------------
    print("Numerical stability")
    amplitudes = []
    for dt in dt_values:
        params.dt = dt
        history, _ = cached_simulation(params, base_method)
        A = max_amplitude_from_history(history, base_method)
        amplitudes.append(A)

    amplitudes = np.array(amplitudes)
    results["stability"]["dt"] = np.array(dt_values)
    results["stability"]["amplitudes"] = amplitudes
    plot_stability(dt_values, amplitudes, base_method)
    plt.show()
    # ---------------------------------------------------------
    # 4. Runtime vs dt
    # ---------------------------------------------------------
    print("Runtime")
    runtimes = []
    for dt in dt_values:
        params.dt = dt
        # force recomputation for runtime measurement (ignore cache)
        history, runtime = run_simulation(params)
        _SIM_CACHE[params.as_key(base_method)] = (history, runtime)
        runtimes.append(runtime)

    runtimes = np.array(runtimes)
    results["runtime"]["dt"] = np.array(dt_values)
    results["runtime"]["runtime"] = runtimes
    plot_runtime(dt_values, runtimes)
    plt.show()
    return results

if __name__ == "__main__":
    for i in ['rk4', 'CN', 'Verlet']:
        run_error_analysis(i)