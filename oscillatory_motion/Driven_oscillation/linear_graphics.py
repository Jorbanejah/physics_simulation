import matplotlib.pyplot as plt
import numpy as np
from Driven_oscillation import DrivenOscillation  
from dataclasses import dataclass
from time import time

##
# ------- Instances -----
##
@dataclass
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

##
# ------- Graphics -----
##

def plot_regime_summary(history, analytical, F0_val, gamma_val):

    """
    Oscillator motion (3x3): 
    Row 1: trajectory (t,q), Row 2: phase diagram (q,v), Row 3: energies
    Columns: (RK4, CN, Verlet) + Analytical
    """

    fig, axes = plt.subplots(4, 3, figsize=(15, 12), tight_layout = True)
    methods = ['rk4', 'CN', 'Verlet']
    
    # Colors for methods
    colors = {'rk4': 'blue', 'CN': 'red', 'Verlet': 'green'}
    
    for i, method in enumerate(methods):
        # Row 1: Position vs time
        axes[0, i].plot(history[method]['t'], history[method]['q'], color=colors[method], linewidth=2, label=method, alpha=0.8)
        axes[0, i].plot(analytical['t'], analytical['x'], 'k--', linewidth=2, label='Analytical', alpha=0.8)

        axes[0, i].set_title('Position')
        axes[0, i].set_xlabel('Time (s)')
        axes[0, i].set_ylabel('Position (mass)')
        axes[0, i].legend()
        axes[0, i].grid(True, alpha=0.3)

        # Row 2: Phase portrait
        axes[1, i].plot(history[method]['q'], history[method]['v'], color=colors[method], linewidth=2, label=method)
        axes[1, i].plot(analytical['x'], analytical['v'], 'k--', linewidth = 2, label = 'Analytical', alpha = 0.8)
        axes[1, i].set_title('Phase Portrait')
        axes[1, i].set_xlabel('Position ')
        axes[1, i].set_ylabel('Velocity ')
        axes[1, i].legend()
        axes[1, i].grid(True, alpha=0.3)

        # Row 3-4: Energies

        E_total = np.array(history[method]['Ek']) +  np.array(history[method]['Ep']) + np.array(history[method]['Wp_diss']) - np.array(history[method]['Wp_drive'])

        axes[2, i].plot(history[method]['t'], E_total)
        axes[2, i].plot(history[method]['t'], history[method]['Ek'], color='red', alpha=0.7, label='Kinetic')
        axes[2, i].plot(history[method]['t'], history[method]['Ep'], color='orange', alpha=0.7, label='Potential')
        
        axes[2, i].set_title ('Energies')
        axes[2, i].set_xlabel('Time (s)')
        axes[2, i].set_ylabel('Energy (J)')
        axes[2, i].legend()
        axes[2, i].grid(True, alpha=0.3)


        axes[3, i].plot(history[method]['t'], history[method]['Wp_diss'], color='red', alpha=0.7, label='Dissipative work')
        axes[3, i].plot(history[method]['t'], history[method]['Wp_drive'], color='orange', alpha=0.7, label='Drive work')
        
        axes[3, i].set_title ('Energies')
        axes[3, i].set_xlabel('Time (s)')
        axes[3, i].set_ylabel('Energy (J)')
        axes[3, i].legend()
        axes[3, i].grid(True, alpha=0.3)
    #plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\regime_summary{F0_val}_{gamma_val}.png", dpi=300, bbox_inches='tight')
# Analitical solution
def beta_vs_amplitude(params, omega0):

    """Study: amplitude vs damping parameter beta for different driving frequencies"""

    betas = np.array([0.1*omega0, 0.2*omega0, 0.3*omega0, 0.5*omega0, 0.8*omega0])
    
    omegas = np.linspace(0.1, 3 * omega0, 50)

    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(betas)) for i in range(len(betas))]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    for i, beta in enumerate(betas):
        amplitude = []

        for omega in omegas:
            denom = np.sqrt((omega0**2 - omega**2)**2 + (2*beta*omega)**2)
            amplitude.append(params.F0/denom)

        label = rf"$\beta$ = {beta/omega0:.1f} $\omega_0$"
        ax1.plot(omegas/omega0, amplitude, color=colors[i], lw=2.5, label=label)
        
    ax1.axvline(1, color='k', ls='--', lw=1, alpha=0.7, label=r"$\omega_0$")
    ax1.axhline(params.F0/omega0**2, color = 'k', ls = '--', lw = 1, alpha = 0.7, label = rf'$F0 / \omega_0 ^2$')
    ax1.set_xlabel(r"$\omega / \omega_0$")
    ax1.set_ylabel("A")
    ax1.set_title("Amplitude vs external frequency _ steady-state")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    max_amplitude = []
    omega_res = []  
    
    for beta in betas:
        # Resonance: 
        omega_r = np.sqrt(max(omega0**2 - 2*beta**2, 0.1))
        omega_res.append(omega_r)
        
        # Maximun power in resonance
        denom_res = np.sqrt((omega0**2 - omega_r**2)**2 + (2*beta*omega_r)**2)
        A_max = params.F0 / denom_res
        max_amplitude.append(A_max)
    
    ax2.semilogy(np.array(betas)/omega0, max_amplitude, 'ro-', markersize=8, lw=2)
    ax2.set_xlabel(r"$\beta / \omega_0$")
    ax2.set_ylabel("A")
    ax2.set_title("Maximum power in resonance vs beta")
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\beta_vs_amplitude.png", dpi=300, bbox_inches='tight')

def beta_vs_power(params, omega0):
    """Study average dissipated power vs external frequency omega, varying beta"""
    
    betas = np.array([0.1*omega0, 0.2*omega0, 0.3*omega0, 0.5*omega0, 0.8*omega0])
    
    omegas = np.linspace(0.1, 3*omega0, 50)
    
    # Colormap
    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(betas)) for i in range(len(betas))]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 1. Main graphic: <P> vs omega for different damping parameter
    for i, beta in enumerate(betas):
        avg_powers = []
        
        for omega in omegas:
            # Theorical solution for average power in steady-state

            denom = (omega0**2 - omega**2)**2 + (2*beta*omega)**2
            P_avg = 0.5 * (params.F0/ params.mass)**2 * (2*beta*omega**2) / denom  
            
            avg_powers.append(P_avg)
        
        label = rf"$\beta$ = {beta/omega0:.1f} $\omega_0$"
        ax1.plot(omegas/omega0, avg_powers, color=colors[i], lw=2.5, label=label)
    
    ax1.axvline(1, color='k', ls='--', lw=1, alpha=0.7, label=r"$\omega_0$")
    ax1.set_xlabel(r"$\omega / \omega_0$")
    ax1.set_ylabel(r"$\langle P \rangle$ (W)")
    ax1.set_title("Average dissipative power vs external frequency")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Secondary Graphics: Maximun power vs beta
    max_powers = []
    omega_res = []  
    
    for beta in betas:
        # Resonance: ω_res ≈ sqrt(ω₀² - 2β²)
        omega_r = np.sqrt(max(omega0**2 - 2*beta**2, 0.1))
        omega_res.append(omega_r)
        
        # Maximun power in resonance
        denom_res = (omega0**2 - omega_r**2)**2 + (2*beta*omega_r)**2
        P_max = 0.5 * (params.F0/params.mass)**2 * (2*beta*omega_r**2) / denom_res
        max_powers.append(P_max)
    
    ax2.semilogy(np.array(betas)/omega0, max_powers, 'ro-', markersize=8, lw=2)
    ax2.set_xlabel(r"$\beta / \omega_0$")
    ax2.set_ylabel(r"$P_{max}$ (W)")
    ax2.set_title("Maximum power in resonance vs beta")
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\beta_vs_power.png", dpi=300, bbox_inches='tight')
    
    return betas, omegas, max_powers

def phase_and_quality_factor():
    
    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k / params.mass) 
    
    betas = np.array([0.05*omega0, 0.1*omega0, 0.2*omega0, 0.4*omega0, 0.7*omega0])
    omegas = np.linspace(0.01, 4*omega0, 200)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(betas)) for i in range(len(betas))]
    
    # 1. Phase
    for i, beta in enumerate(betas):
        phases = []
        for omega in omegas:

            phase = np.arctan2(2*beta*omega, omega0**2 - omega**2)
            phases.append(np.degrees(phase))  # Degrees
    
        Q = omega0 / (2*beta)  # Quality factor
        label = rf"Q={Q:.1f} ($\beta$={beta/omega0:.2f}$\omega_0$)"
        axes[0,0].plot(omegas/omega0, phases, color=colors[i], lw=3, label=label)
    
    axes[0,0].axvline(1.0, color='k', ls='--', lw=2, alpha=0.8, label=r"$\omega = \omega_0$")
    axes[0,0].axhline(90, color='blue', ls=':', lw=2, alpha=0.7, label=r"$\phi$=90°")
    axes[0,0].set_xlabel(r"$\omega/\omega_0$")
    axes[0,0].set_ylabel(r"$\phi$ (grados)")
    axes[0,0].set_title('Phase through resonance')
    axes[0,0].legend(framealpha=0.95)
    axes[0,0].grid(True, alpha=0.3)
    
    # 2. Quality factor vs width
    delta_omegas = []
    for beta in betas:

        Q = omega0 / (2*beta)

        delta_omega = omega0 / Q
        delta_omegas.append(delta_omega)
    
    axes[0,1].semilogy(np.array(betas)/omega0, [omega0/(2*b) for b in betas], 'ro-', markersize=10, lw=3, label='Q=ω₀/(2β)')
    axes[0,1].semilogy(np.array(betas)/omega0, delta_omegas/omega0, 'bs-', markersize=8, lw=3, label='Δω/ω₀ = 1/Q')
    axes[0,1].set_xlabel(r"$\beta/\omega_0$")
    axes[0,1].set_ylabel("Adimensional factor")
    axes[0,1].set_title(r"Quality factor vs width $\delta \omega$")
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # 3. φ vs ω normalization for Q
    for i, beta in enumerate(betas):
        Q = omega0 / (2*beta)
        phi_Q = []
        for omega in omegas:
            phase = np.arctan2(2*beta*omega, omega0**2 - omega**2)
            phi_Q.append(np.degrees(phase))
        
        axes[1,0].plot(omegas/omega0, phi_Q, color=colors[i], lw=2.5)
    
    axes[1,0].set_xlabel(r"$\omega / \omega_0$")
    axes[1,0].set_ylabel(r"$\phi$ (grados)")
    axes[1,0].set_title("Phase vs omega normalize for Q")
    axes[1,0].grid(True, alpha=0.3)
    
    # 4. Curva maestra universal φ vs ω/ω₀ (independiente de Q)
    u = np.linspace(0.1, 5, 200)  # Represent ω/ω₀
    phi_universal = np.degrees(np.arctan2(2*u, 1 - u**2)) 
    axes[1,1].plot(u, phi_universal, 'k-', lw=4, label="Q→∞")
    axes[1,1].set_xlabel(r"$\omega/\omega_0$")
    axes[1,1].set_ylabel(r"$\phi$ (grados)")
    axes[1,1].set_title("Limit  Q→∞")
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\phase_and_qualityfactor.png", dpi=300, bbox_inches='tight')
    
    # Results table
    print("\n" + "="*60)
    print("Quality factor:")
    print("="*60)
    for beta in betas:
        Q = omega0 / (2*beta)
        delta_w = omega0 / Q
        print(f"β/ω₀={beta/omega0:.2f} → Q={Q/omega0:.2f} → Δω/ω₀={delta_w/omega0:.3f}")
    print("="*60)

##  
# ------- Main execution -------
##
if __name__ == "__main__":

    # Create linear oscillator
    params = DrivenOscillationParams()
    
    numerical_method = ['rk4', 'CN', 'Verlet']
    omega0 = np.sqrt(params.k/params.mass)
    
    F0_list = [0, 0, 2.0]
    gamma_list = [0, 1, 2]

    #Remenber: w_d = np.sqrt(k/m - gamma/(2m)) and that must be real

    for idx, (F0_val, gamma_val) in enumerate(zip(F0_list, gamma_list)):

        print(f"Running F0={F0_val}, gamma={gamma_val}")
    
        osc = DrivenOscillation( q0=params.q0, dq0=params.dq0, m=params.mass, gamma=gamma_val, F0=F0_val, omega=params.omega, t_max=params.t_max, dt=params.dt, system='linear', k=params.k, F_external=params.F_external)
    
        model = osc.run() 
    
        plot_regime_summary(model.history, model.analytical, F0_val, gamma_val) 
    
    # Resonance curve
    #beta_vs_power(params, omega0)
    #beta_vs_amplitude(params, omega0)

    #Quality factor

    #phase_and_quality_factor()

    plt.show()