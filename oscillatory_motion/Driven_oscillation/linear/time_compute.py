import numpy as np
import time
import matplotlib.pyplot as plt
from Driven_oscillation import DrivenOscillation

class DrivenOscillationParams:
    """Parameters for linear driven oscillator"""
    mass: float = 2.0
    gamma: float = 1.0
    k: float = 1.0

    omega: float = 2.0
    F0: float = 2.0
    F_external: str = 'cos'

    q0: float = 2.0
    dq0: float = 2.0

    dt: float = 0.01
    t_max: float = 25.0
    
    system: str = 'linear'

def energy_surface_omega(params, omega0, methods, dts):

    betas = np.linspace(0*omega0, 1.5*omega0, 10)
    omegas = np.linspace(0.1 * omega0, 3 * omega0, 10)


    results = {
        method: {
            "dt": [],
            "betas": betas,
            "omegas": omegas,
            "drift_maps": [],
            "compute_times": []
        }
        for method in methods
    }

    for method in methods:
        print(f"\n=== Running method: {method} ===")

        for dt in dts:
            print(f"  → dt = {dt}")

            drift_map = np.zeros((len(betas), len(omegas)))

            t0 = time.perf_counter()

            for i, beta in enumerate(betas):
                for j, omega in enumerate(omegas):

                    cycles = 30
                    T = 2*np.pi / abs(omega)
                    t_max = cycles * T

                    osc = DrivenOscillation(
                        q0=params.q0, dq0=params.dq0, m=params.mass,
                        gamma=beta * 2 * params.mass,
                        F0=params.F0, omega=omega, t_max=t_max, dt=dt,
                        system='linear', k=params.k, F_external=params.F_external
                    )

                    model = osc.run()
                    history = model.history

                    E_total = (
                        np.array(history[method]['Ek']) +
                        np.array(history[method]['Ep']) +
                        np.array(history[method]['Wp_diss']) -
                        np.array(history[method]['Wp_drive'])
                    )

                    drift_map[i, j] = (np.max(E_total) - np.min(E_total)) / np.mean(E_total)

            t1 = time.perf_counter()
            compute_time = t1 - t0

            # store results
            results[method]["dt"].append(dt)
            results[method]["drift_maps"].append(drift_map)
            results[method]["compute_times"].append(compute_time)

    return results


def plot_compute_times(results):
    plt.figure(figsize=(8,5))

    for method, data in results.items():
        plt.plot(data["dt"], data["compute_times"], marker='o', label=method)

    plt.xlabel("dt")
    plt.ylabel("Compute time (s)")
    plt.title("Compute time vs dt")
    plt.legend()
    plt.grid(True)
    plt.savefig("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\time_vs_dt.png", dpi=300, bbox_inches='tight')

def plot_drift_vs_dt(results):
    plt.figure(figsize=(8,5))

    for method, data in results.items():
        avg_drift = [np.mean(drift_map) for drift_map in data["drift_maps"]]
        plt.plot(data["dt"], avg_drift, marker='o', label=method)

    plt.xlabel("dt")
    plt.ylabel("Average drift")
    plt.title("Energy drift vs dt")
    plt.legend()
    plt.grid(True)
    plt.savefig("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\drift_energy_vs_dt.png", dpi=300, bbox_inches='tight')

def plot_drift_surface(results, method, index_dt):
    betas = results[method]["betas"]
    omegas = results[method]["omegas"]
    drift_map = results[method]["drift_maps"][index_dt]

    B, W = np.meshgrid(betas, omegas, indexing='ij')

    fig = plt.figure(figsize=(12,5))

    ax = fig.add_subplot(121, projection='3d')
    ax.plot_surface(B, W, drift_map, cmap='viridis') 
    ax.set_xlabel(r"$\beta / \omega_0$")
    ax.set_ylabel(r"$\omega / \omega_0$")
    ax.set_zlabel(r"$\Delta E_{total}$")
    ax.set_title(f"{method} drift surface (dt={results[method]['dt'][index_dt]})")

    ax2 = fig.add_subplot(122)
    im = ax2.imshow(drift_map, origin='lower', cmap='viridis', extent=[omegas[0] / omega0, omegas[-1] / omega0, betas[0] / omega0, betas[-1] / omega0], aspect = 'auto') #Normalize with omega0 the drift map
    ax2.set_xlabel(r"$\omega / \omega_0$")
    ax2.set_ylabel(r"$\beta / \omega_0$")
    plt.colorbar(im, ax=ax2)
    ax2.set_title("Heatmap")

    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\heatmap_errors{method}_{index_dt}.png", dpi=300, bbox_inches='tight')

if __name__ == "__main__":

    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k/params.mass)

    methods = ["rk4", "CN", "Verlet"]
    dts = [0.01, 0.05, 0.1, 0.5]

    results = energy_surface_omega(params, omega0, methods, dts)

    #plot_compute_times(results)

    #plot_drift_vs_dt(results)
    
# index_dt = 0 - 0.01, 1 - 0.05, 2 - 0.1, 3 - 0.5
    for i in methods:
        plot_drift_surface(results, method= i, index_dt = 0)
        plot_drift_surface(results, method= i, index_dt= 3)

    plt.show()