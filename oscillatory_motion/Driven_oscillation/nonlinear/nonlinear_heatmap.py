import numpy as np
import matplotlib.pyplot as plt
from Driven_oscillation import DrivenOscillation
from dataclasses import dataclass

@dataclass
class DrivenOscillationParams(): 
        #Pendulum
        mass: float = 2
        gamma: float = 1
        L: float = 2

        # External forces
        omega: float = 2
        F0: float = 1
        F_external: str = 'cos'
        system: str = 'nonlinear'

        #Initial condition
        q0: float = np.deg2rad(30)
        dq0: float = 0 

        #Times
        dt: float = 0.01
        t_max: float = 20

def energy_surface_omega(params, omega0, method='rk4'):

    # Parameter sweeps
    betas = np.linspace(0*omega0, 1.5*omega0, 40)  # Physical damping
    omegas = np.linspace(0.1 * omega0, 3 * omega0, 40)
    dt = [0.01, 0.1, 0.5]
    
    for dt in dt:
        # Storage for drift values
        drift_map = np.zeros((len(betas), len(omegas)))

        for i, beta in enumerate(betas):
            for j, omega in enumerate(omegas):
                #We want that every omega contribute equally
                cycles = 30
                T = 2*np.pi / abs(omega)
                t_max = cycles * T

                osc = DrivenOscillation(q0=params.q0, dq0=params.dq0, m=params.mass, gamma= beta * params.L**2 * params.mass, F0=params.F0, omega=omega, t=t_max, dt=dt, system='nonlinear', L=params.L, F_external=params.F_external)
                model = osc.run()
                history = model.history

               
                E_total = np.array(history[method]['Ek']) +  np.array(history[method]['Ep']) + np.array(history[method]['Wp_diss']) - np.array(history[method]['Wp_drive'])

                dE = (np.max(E_total) - np.min(E_total)) / np.mean(E_total)

                drift_map[i, j] = dE   # scalar drift measure

        print(f'Running {method=} with {dt=}')

        # ---- 3D SURFACE PLOT ----
        B, W = np.meshgrid(betas / omega0, omegas / omega0, indexing='ij')

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(121, projection='3d')

        surf = ax.plot_surface(B, W, drift_map, cmap='viridis', edgecolor='none')
        ax.set_xlabel(r'$\beta / \omega_0$')
        ax.set_ylabel(r'$\omega / \omega_0$')
        ax.set_zlabel(r'$\max |E_{total}|$')
        ax.set_title(f'Energy Drift Surface ({method})')

        # ---- 2D HEATMAP ----
        ax2 = fig.add_subplot(122)
        im = ax2.imshow(drift_map, extent=[omegas[0]/omega0, omegas[-1]/omega0, betas[0]/omega0, betas[-1]/omega0], origin='lower', aspect='auto', cmap='viridis')
        ax2.set_xlabel(r'$\omega / \omega_0$')
        ax2.set_ylabel(r'$\beta / \omega_0$')
        ax2.set_title(f'Energy Drift Map ({method})')
        fig.colorbar(im, ax=ax2)

        plt.tight_layout()
        plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\nonlinear_heatmap_{method}_{dt}.png", dpi=300, bbox_inches='tight')
    
def energy_surface_F0(params, omega0, method='rk4'):

    # Parameter sweeps
    betas = np.linspace(0*omega0, 1.5*omega0, 20)  # Physical damping
    F0s = np.linspace(0.1 * omega0, 3 * omega0, 20)
    dt = [0.01, 0.05, 0.1, 0.5]
    
    for dt in dt:
        # Storage for drift values
        drift_map = np.zeros((len(betas), len(F0s)))

        for i, beta in enumerate(betas):
            for j, F0 in enumerate(F0s):
                #We want that every omega contribute equally
                cycles = 30
                T = 2*np.pi / abs(F0)
                t_max = cycles * T

                osc = DrivenOscillation(q0=params.q0, dq0=params.dq0, m=params.mass,gamma=beta * 2 * params.mass, F0=F0, omega=params.omega, t_max=t_max, dt=dt, system='nonlinear', k=params.k, F_external=params.F_external)
                model = osc.run()
                history = model.history

                E_total = np.array(history[method]['Ek']) +  np.array(history[method]['Ep']) + np.array(history[method]['Wp_diss']) - np.array(history[method]['Wp_drive'])

                dE = (np.max(E_total) - np.min(E_total)) / np.mean(E_total)

                drift_map[i, j] = dE   # scalar drift measure

        print(f'Running {method=} with {dt=}')

        # ---- 3D SURFACE PLOT ----
        B, W = np.meshgrid(betas / omega0, F0s / omega0, indexing='ij')

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(121, projection='3d')

        surf = ax.plot_surface(B, W, drift_map, cmap='viridis', edgecolor='none')
        ax.set_xlabel(r'$\beta / \omega_0$')
        ax.set_ylabel(r'$F0 / \omega_0$')
        ax.set_zlabel(r'$\max |E_{total}|$')
        ax.set_title(f'Energy Drift Surface ({method})')

        # ---- 2D HEATMAP ----
        ax2 = fig.add_subplot(122)
        im = ax2.imshow(drift_map, extent=[F0s[0]/omega0, F0s[-1]/omega0, betas[0]/omega0, betas[-1]/omega0], origin='lower', aspect='auto', cmap='viridis')
        ax2.set_xlabel(r'F0 / \omega_0$')
        ax2.set_ylabel(r'$\beta / \omega_0$')
        ax2.set_title(f'Energy Drift Map ({method})')
        fig.colorbar(im, ax=ax2)

        plt.tight_layout()
        #plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\energy_surface_F0_{method}_{dt}.png", dpi=300, bbox_inches='tight')
    
if __name__ == "__main__":

    params = DrivenOscillationParams()
    g = 9.81
    omega0 = np.sqrt(g/params.L)

    energy_surface_omega(params, omega0)
    plt.show()