from dataclasses import dataclass
import numpy as np
from Driven_oscillation import DrivenOscillation
'''
The following code store the differents trajectories that will be used inside the animations: Poincare sections, Lypunov coefficient and bifurcation diagram.
We have used the Taylor parameters, and stored the trajectories inside of trajectories.npz:


Taylor forms are different too: 

d2(theta) + beta*d(theta) + sin(theta)  = gamma * cos(omega)

𝛽 = B/(𝑚 * L**2 * g/L)
gamma = F0 /(𝑚 * g * L)
omega = Ω/omega0 ---> Ω omega driven


1.060 < gamma < 1.087

omega_driven = 2/3

beta = 0.25

q0 = 0
dq0 =0 
omega_0 = 1 --- omega_0 = sqrt(g/L) then L = 9.8

Taylor's book does not explicity take the steady state, however, he discards 200-300 samples and takes the folowing 100-200 samples
So, period(T) = 2pi/omega = 3pi -----> t = ?

'''
@dataclass
class DrivenOscillationParameter:
    #Initial condition

    q0: float = 0
    dq0: float = 0

    #
    q02: float = np.deg2rad(0.1)

    #Innate parameter
    L: float = 1
    m: float = 2
    gamma: float = 0.1

    #External force
    F0: float = 2
    F_external: str = 'cos'
    omega: float = 2 * np.pi * 0.5

    #Time
    t: int = 800 #We need a huge time even 2000
    dt: float = 0.01

    #System
    system: str = 'nonlinear'

class DynamicalSystem:
    """
    Encapsulates numerical simulation for a sweep of the forcing amplitude.
    """

    def __init__(self):
        self.params = DrivenOscillationParameter()
        self.alphas = np.linspace(1.060, 1.087, 40)  # Range showing normal→chaotic
        self.methods = ["rk4"]

        # Allocate storage
        self.poincare_sections = {
            method: {"q": {}, "dq": {}, "q_full": {}, "dq_full": {}} for method in self.methods
        }

    def extract_poincare_section(self, q_traj, dq_traj, omega, t_start=400):
        """Extract Poincaré section by sampling at driving period."""

        period = 2 * np.pi / omega
        dt = self.params.dt
        steps_per_period = int(period/dt)

        t_sample = np.arange(int(t_start/dt), len(q_traj), steps_per_period)
        
        return q_traj[t_sample], dq_traj[t_sample]

    def run_parameter_sweep(self):
        """Compute Poincaré sections for all alpha values."""
        print("Computing Poincaré sections...")
        
        for i, alpha in enumerate(self.alphas):
            progress = (i + 1) / len(self.alphas)
            bar_length = 12
            filled = int(progress * bar_length)
            bar = "█" * filled + "-" * (bar_length - filled)

            print(rf"[{bar}]  {progress*100:5.1f}%   $\alpha$ = {alpha:.2f}", end="\r", flush=True)
            
            # Update F0 = alpha
            F0 = alpha * self.params.m * self.params.L**2
            
            # Create oscillator with current parameters
            osc = DrivenOscillation(q0=self.params.q0, dq0=self.params.dq0, m=self.params.m, gamma=self.params.gamma, F0=F0, omega=self.params.omega, t = self.params.t, system=self.params.system, L=self.params.L, F_external = self.params.F_external)
            
            # Run simulation
            model = osc.run()
            
            # Extract Poincaré sections for each method
            for method in self.methods:
                q_full = np.array(model.history[method]["q"])
                dq_full = np.array(model.history[method]["v"])

                # Extract Poincaré section (discard transients)
                q_poincare, dq_poincare = self.extract_poincare_section(q_full, dq_full, self.params.omega)
                
                self.poincare_sections[method]["q"][alpha] = q_poincare
                self.poincare_sections[method]["dq"][alpha] = dq_poincare
                self.poincare_sections[method]["q_full"][alpha] = q_full
                self.poincare_sections[method]["dq_full"][alpha] = dq_full
        return self.poincare_sections

    def store(self, filename="C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz"):
        """Store results."""
        
        np.savez(filename, **self.poincare_sections)
        print(f"Data stored in {filename}")

# Pre-compute data (run this first)
if __name__ == "__main__":
    system = DynamicalSystem()
    sections = system.run_parameter_sweep()
    system.store()
    data = np.load("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz", allow_pickle=True)
    rk4 = data["rk4"].item()

    for alpha, arr in rk4["q"].items():
        print(alpha, len(arr), arr)
