
import matplotlib.pyplot as plt
import numpy as np
from Driven_oscillation import DrivenOscillation  
from dataclasses import dataclass
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
    F0: float = 6.0
    F_external: str = 'cos'

    q0: float = 2.0
    dq0: float = 0.0

    dt: float = 0.01
    t_max: float = 25.0
    
    system: str = 'linear'

##
# ------- Graphics -----
##

def plot_regime_summary(history, analytical):

    """
    Oscillator motion (3x3): 
    Row 1: trajectory (t,q), Row 2: phase diagram (q,v), Row 3: energies
    Columns: RK4, CN, Verlet + Analytical
    """

    fig, axes = plt.subplots(4, 3, figsize=(15, 10), tight_layout = True)
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
        axes[1, i].set_title('Phase Portrait')
        axes[1, i].set_xlabel('Position (mass)')
        axes[1, i].set_ylabel('Velocity (mass/s)')
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
    plt.show()

# Analitical solution
def beta_vs_amplitude(betas, omegas):
    """Study amplitude vs damping parameter beta for different driving frequencies"""
    params = DrivenOscillationParams()
    
    fig, axes = plt.subplots(1, len(omegas), figsize=(5*len(omegas), 4))
    if len(omegas) == 1:
        axes = [axes]
    
    for idx, (beta, omega) in enumerate(zip(betas, omegas)):
        ax = axes[idx] if len(omegas) > 1 else axes
        
        amplitudes = []
        for b in betas:
            # Update parameters
            temp_params = DrivenOscillationParams(
                gamma=2*params.mass*b, omega=omega, **{k: v for k, v in params.__dict__.items() 
                                                      if k != 'gamma' and k != 'omega'}
            )
            
            osc = DrivenOscillation(
                q0=temp_params.q0, dq0=temp_params.dq0, mass=temp_params.mass,
                gamma=temp_params.gamma, F0=temp_params.F0, omega=temp_params.omega,
                t=temp_params.t_max, dt=temp_params.dt, system='linear',
                k=params.k, F_external=temp_params.F_external
            )
            
            model = osc.run()
            history, _ = model.run()
            
            # Compute steady-state amplitude (last 20% of trajectory)
            n_steady = int(0.8 * len(history['rk4']['q']))
            amp = np.max(np.abs(history['rk4']['q'][n_steady:]))
            amplitudes.append(amp)
        
        ax.plot(betas, amplitudes, 'o-', linewidth=2, markersize=6)
        ax.set_xlabel('Damping coefficient β')
        ax.set_ylabel('Steady-state amplitude')
        ax.set_title(f'ω = {omega:.1f} rad/s')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def beta_vs_velocity(betas, omega=2.0):
    """Study maximum velocity vs damping parameter beta"""
    params = DrivenOscillationParams()
    params.omega = omega
    
    fig, ax = plt.subplots(figsize=(8, 6))
    max_velocities = []
    
    for beta in betas:
        temp_params = DrivenOscillationParams(gamma=2*params.mass*beta, omega=omega)
        
        osc = DrivenOscillation(
            q0=temp_params.q0, dq0=temp_params.dq0, mass=temp_params.mass,
            gamma=temp_params.gamma, F0=temp_params.F0, omega=temp_params.omega,
            t=temp_params.t_max, dt=temp_params.dt, system='linear',
            k=params.k, F_external=temp_params.F_external
        )
        
        model = osc.run()
        history, _ = model.run()
        
        # Steady-state max velocity
        n_steady = int(0.8 * len(history['rk4']['v']))
        max_v = np.max(np.abs(history['rk4']['v'][n_steady:]))
        max_velocities.append(max_v)
    
    ax.plot(betas, max_velocities, 's-', color='green', linewidth=2, markersize=8)
    ax.set_xlabel('Damping coefficient β')
    ax.set_ylabel('Maximum steady-state velocity (m/s)')
    ax.set_title(f'Maximum velocity vs β (ω = {omega} rad/s)')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def beta_vs_power():
    """Study average dissipated power vs external frequency ω, varying β"""
    
    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k / params.mass)  # Frecuencia natural CORRECTA
    
    # Rango de β (fracciones de ω₀)
    betas = np.array([0.1*omega0, 0.2*omega0, 0.3*omega0, 0.5*omega0, 0.8*omega0])
    
    # Rango de frecuencias externas
    omegas = np.linspace(0.1, 3*omega0, 50)
    
    # Colormap
    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(betas)) for i in range(len(betas))]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 1. GRÁFICA PRINCIPAL: <P> vs ω para diferentes β
    for i, beta in enumerate(betas):
        avg_powers = []
        
        for omega in omegas:
            # FÓRMULA ANALÍTICA EXACTA para potencia promedio steady-state
            denom = (omega0**2 - omega**2)**2 + (2*beta*omega)**2
            P_avg = 0.5 * (params.F0/ params.mass)**2 * (2*beta*omega**2) / denom  # ← FÓRMULA CORRECTA
            
            avg_powers.append(P_avg)
        
        label = rf"$\beta$ = {beta/omega0:.1f} $\omega_0$"
        ax1.plot(omegas/omega0, avg_powers, color=colors[i], lw=2.5, label=label)
    
    ax1.axvline(1.0, color='k', ls='--', lw=1, alpha=0.7, label=r"$\omega_0$")
    ax1.set_xlabel(r"$\omega / \omega_0$")
    ax1.set_ylabel(r"$\langle P \rangle$ (W)")
    ax1.set_title("Potencia disipada promedio vs frecuencia externa")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. GRÁFICA SECUNDARIA: Máxima potencia vs β
    max_powers = []
    omega_res = []  # Frecuencia de resonancia
    
    for beta in betas:
        # Resonancia: ω_res ≈ sqrt(ω₀² - 2β²)
        omega_r = np.sqrt(max(omega0**2 - 2*beta**2, 0.1))
        omega_res.append(omega_r)
        
        # Potencia máxima en resonancia
        denom_res = (omega0**2 - omega_r**2)**2 + (2*beta*omega_r)**2
        P_max = 0.5 * (params.F0/params.mass)**2 * (2*beta*omega_r**2) / denom_res
        max_powers.append(P_max)
    
    ax2.semilogy(np.array(betas)/omega0, max_powers, 'ro-', markersize=8, lw=2)
    ax2.set_xlabel(r"$\beta / \omega_0$")
    ax2.set_ylabel(r"$P_{max}$ (W)")
    ax2.set_title("Potencia máxima de resonancia vs β")
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return betas, omegas, max_powers

def beta_vs_power_numerical():
    """Validación NUMÉRICA vs ANALÍTICA"""
    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k / params.mass)
    betas = np.array([0.2*omega0, 0.4*omega0])
    omega_test = omega0  # Frecuencia natural
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    for i, beta in enumerate(betas):
        ax = axes[i]
        
        # Simulación numérica
        temp_params = DrivenOscillationParams(gamma=2*params.mass*beta)
        osc = DrivenOscillation(q0= temp_params.q0, dq0=temp_params.dq0, m=temp_params.mass, gamma=temp_params.gamma, F0=temp_params.F0, omega=omega_test, t=50, dt=0.01, system='linear', k=params.k, F_external='sin')
        model = osc.run()
        
        t = np.array(model.history['rk4']['t'])
        P_drive = np.array(model.history['rk4']['Wp_drive'])
        P_diss = np.array(model.history['rk4']['Wp_diss'])
        
        # Potencia promedio steady-state (últimos 20%)
        n_steady = int(0.8 * len(t))
        P_avg_num = (P_diss[-1] - P_diss[n_steady]) / (t[-1] - t[n_steady])
        
        # Analítica
        denom = (omega0**2 - omega_test**2)**2 + (2*beta*omega_test)**2
        P_avg_ana = 0.5 * (params.F0/params.mass)**2 * (2*beta*omega_test**2) / denom
        
        print(f"β={beta/omega0:.1f}ω₀: P_num={P_avg_num:.4f}, P_ana={P_avg_ana:.4f}, Error={100*abs(P_avg_num-P_avg_ana)/P_avg_ana:.1f}%")
        
        ax.plot(t[n_steady:], P_drive[n_steady:], 'g-', alpha=0.7, label='P_drive')
        ax.plot(t[n_steady:], P_diss[n_steady:], 'r-', label='P_diss')
        ax.axhline(P_avg_num, color='k', ls='--', label=f'P_avg={P_avg_num:.3f}')
        ax.set_title(rf'β = {beta/omega0:.1f}ω₀ (ω={omega_test/omega0:.1f}ω₀)')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()


##
# ------- Main execution -------
##
if __name__ == "__main__":

    # Create linear oscillator
    params = DrivenOscillationParams()
    osc = DrivenOscillation(q0=params.q0, dq0=params.dq0, m=params.mass, gamma=params.gamma, F0=params.F0, omega=params.omega, t=params.t_max, dt=params.dt, system='linear',  k = params.k, F_external=params.F_external)
    
    model = osc.run()
    history, analytical = model.run()
    
    # Plot regime summary
    #plot_regime_summary(history, analytical)
    
    # Parametric studies
    betas = np.linspace(0.1, 1.0, 100)
    #omegas = [1.0, 1.5, 2.0]
    
    #beta_vs_amplitude(betas, omegas)
    #beta_vs_velocity(betas)
    beta_vs_power()
    beta_vs_power_numerical()

# Poincare section (2-d and 3-d)

# Bifurcation diagrams

