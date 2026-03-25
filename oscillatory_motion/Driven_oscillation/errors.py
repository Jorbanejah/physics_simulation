import numpy as np
import matplotlib.pyplot as plt
from Driven_oscillation import DrivenOscillation
from dataclasses import dataclass

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

def beta_vs_power_numerical():
    """Numerical vs Analitical validation"""

    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k / params.mass)
    betas = np.array([i* omega0 for _, i in enumerate(np.linspace(0.1 * omega0, 0.8 * omega0, 10))])
    omega_test = omega0  # Natural frequency
    
    numerical_method = ['rk4', 'CN', 'Verlet']
    errors = {name: [] for name in numerical_method}
    plt.figure(figsize= (10, 6))

    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(numerical_method)) for i in range(len(numerical_method))]
    

    for _, beta in enumerate(betas):
        # Analitical
        denom = (omega0**2 - omega_test**2)**2 + (2*beta*omega_test)**2
        P_avg_ana = 0.5 * (params.F0/params.mass)**2 * (2*beta*omega_test**2) / denom

        comparation = {
            name: {'t':[], 'Wp_diss': []} for name in numerical_method 
        }

        for j, methods in enumerate(numerical_method):

            # Numerical simulation
            temp_params = DrivenOscillationParams(gamma=2*params.mass*beta)
            osc = DrivenOscillation(q0= temp_params.q0, dq0=temp_params.dq0, m=temp_params.mass, gamma=temp_params.gamma, F0=temp_params.F0, omega=omega_test, t_rk4=50, dt=0.01, system='linear', k=params.k, F_external='sin')
            model = osc.run()

            comparation[methods]['t'] = np.array(model.history[methods]['t'])
            comparation[methods]['Wp_diss'] = np.array(model.history[methods]['Wp_diss'])
            
            n_steady = int(0.8 * len(comparation[methods]['t']))
            P_avg_num = (comparation[methods]['Wp_diss'][-1] - comparation[methods]['Wp_diss'][n_steady])/( comparation[methods]['t'][-1] - comparation[methods]['t'][n_steady])

            #print(rf"$\beta$={beta/omega0:.1f}$\omega_0$: P_num={P_avg_num:.4f}, P_ana={P_avg_ana:.4f}, Error={100*abs(P_avg_num-P_avg_ana)/P_avg_ana:.1f}%. for {methods}")
            Error = 100*abs(P_avg_num-P_avg_ana)/P_avg_ana
            errors[methods].append(Error)

    for j, method in enumerate(numerical_method):

        plt.plot(betas/omega0, errors[method], '-o', color=colors[j], label=method)

    plt.xlabel(r'$\beta / \omega_0$')
    plt.ylabel('Relative Error (%)')
    plt.title("Absolute Error: numerical vs analitical")
    plt.legend()
    plt.grid(True)
    plt.show()

def beta_vs_amplitude_numerical():

    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k / params.mass)
    betas = np.array([i* omega0 for _, i in enumerate(np.linspace(0.1 * omega0, 0.8 * omega0, 10))])
    
    numerical_method = ['rk4', 'CN', 'Verlet']
    errors = {name: [] for name in numerical_method}

    plt.figure(figsize= (10, 6))

    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(numerical_method)) for i in range(len(numerical_method))]
    omega = omega0 # Resonance

    for _, beta in enumerate(betas):

        denom = np.sqrt((omega0**2 - omega**2)**2 + (2*beta*omega)**2) 
        alpha = params.F0 / params.mass
        amplitude = (alpha/denom)

        comparation = {
            name: {'t':[], 'Wp_diss': []} for name in numerical_method 
        }
        for j, methods in enumerate(numerical_method):
            
            try:
                temp_params = DrivenOscillationParams(gamma=2*params.mass*beta)
                osc = DrivenOscillation(q0= temp_params.q0, dq0=temp_params.dq0, m=temp_params.mass, gamma=temp_params.gamma, F0=temp_params.F0, omega=omega, t_rk4=50, dt=0.01, system='linear', k=params.k, F_external='sin')
                model = osc.run()

                comparation[methods]['t'] = np.array(model.history[methods]['t'])
                comparation[methods]['q'] = np.array(model.history[methods]['q'])
            
                n_steady = int(0.8 * len(comparation[methods]['t']))
                A_avr = (max(comparation[methods]['q'][n_steady:]) - min(comparation[methods]['q'][n_steady:]))

                Error = 100 * abs(A_avr - amplitude) / amplitude
                errors[methods].append(Error)
                #print(rf"$\beta$={beta/omega0:.1f}$\omega_0$: A_num={A_avr:.4f}, A_ana={amplitude:.4f}, Error={Error}%. for {methods}")
            except Exception as e:

                print(f"Error {method}, β={beta/omega0:.2f}: {e}")
                errors[method].append(np.nan)

    for j, method in enumerate(numerical_method):
         
        plt.plot(betas/omega0, errors[method], '-o', color=colors[j], label=method)
    plt.xlabel(r'$\beta / \omega_0$')
    plt.ylabel('Relative Error (%)')
    plt.title("Absolute Error: numerical vs analitical")
    plt.legend()
    plt.grid(True)
    plt.show()

    
def validate_phase_fitting():
    
    params = DrivenOscillationParams()
    omega0 = np.sqrt(params.k / params.mass)
    alpha = params.F0 / params.mass
    
    beta = 0.2 * omega0
    omega_drive = 1.3 * omega0
    
    # Simulation
    osc = DrivenOscillation(
        q0=0.01, dq0=0.0, m=params.mass, gamma=2*params.mass*beta,
        F0=params.F0, omega=omega_drive, t=200, dt=0.002,
        system='linear', k=params.k, F_external='sin'
    )
    model = osc.run()
    
    t = np.array(model.history['rk4']['t'])
    q = np.array(model.history['rk4']['q'])
    n_steady = int(0.95 * len(t))
    
    t_ss, q_ss = t[n_steady:], q[n_steady:]
    F_theory = alpha * np.sin(omega_drive * t_ss)
    
    # Crossing correlation
    from scipy.signal import correlate
    corr = correlate(F_theory, q_ss, mode='full')
    delay = len(F_theory) - np.argmax(corr)
    phi_fit = (delay * 0.002 * omega_drive) * 180 / np.pi  # dt=0.002
    
    # Analitically
    phi_ana = np.degrees(np.arctan2(2*beta*omega_drive, omega0**2 - omega_drive**2))
    
    print(fr"Phase: $\phi_fit$={phi_fit:.1f}°, $\phi_ana$= {phi_ana:.1f}°")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Trajectories
    ax1.plot(t_ss, q_ss, 'b-', lw=2, label='x(t)')
    ax1.plot(t_ss, F_theory, 'r--', lw=2, label='sin(ωt)')
    ax1.legend()
    ax1.set_title('F(t) vs x(t)')
    ax1.grid(True, alpha=0.3)
    
    # Correlation
    corr = correlate(F_theory, q_ss, mode='full')
    ax2.plot(corr, 'g-', lw=2)
    ax2.axvline(len(F_theory)-1, color='r', ls='--', label=f'φ={phi_fit:.1f}°')
    ax2.set_xlabel('Retardo')
    ax2.set_title('Correlación cruzada')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    # Create linear oscillator
    params = DrivenOscillationParams()
    osc = DrivenOscillation(q0=params.q0, dq0=params.dq0, m=params.mass, gamma=params.gamma, F0=params.F0, omega=params.omega, t_rk4=params.t_max, dt=params.dt, system='linear',  k = params.k, F_external=params.F_external)
    
    model = osc.run()
    history, analytical = model.run()
    
    validate_phase_fitting()
    #beta_vs_amplitude_numerical()
    #beta_vs_power_numerical()
    

