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
    L: float = 2.0

    omega: float = 2.0
    F0: float = 1.0
    F_external: str = 'cos'

    theta0: float = np.deg2rad(30)
    omega0: float = 1.0

    dt: float = 0.01
    t_max: float = 20.0
    
    system: str = 'nonlinear'

##
# ------- Graphics -----
##

def plot_regime_summary(history, F0_val, gamma_val):

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

        t = np.array(history[method]['t'])
        Ek = np.array(history[method]['Ek'])
        Ep = np.array(history[method]['Ep'])
        Wp_diss = np.array(history[method]['Wp_diss']) 
        Wp_drive = np.array(history[method]['Wp_drive'])
        
        # Just mechanical energy!
        E_mech = Ek + Ep  
        
        # Net work done ON system
        W_net = Wp_drive - Wp_diss  # Drive adds, dissipation subtracts
        
        # Row 1: Position vs time
        axes[0, i].plot(t, history[method]['q'], color=colors[method], linewidth=2, label=method, alpha=0.8)

        axes[0, i].set_title('Position')
        axes[0, i].set_xlabel('Time (s)')
        axes[0, i].set_ylabel('Position (mass)')
        axes[0, i].legend()
        axes[0, i].grid(True, alpha=0.3)

        # Row 2: Phase portrait
        axes[1, i].plot(history[method]['q'], history[method]['v'], color=colors[method], linewidth=2, label=method)
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
    plt.savefig(f"C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\figures\\regime_summary_nonlinear{F0_val}_{gamma_val}.png", dpi=300, bbox_inches='tight')
def average_power(history, method, discard_fraction=0.75):
    t = np.array(history[method]["t"])
    W = np.array(history[method]["Wp_drive"])

    # discard initial transient
    N = int(len(t) * discard_fraction)
    t_ss = t[N:]
    W_ss = W[N:]

    # average power
    return (W_ss[-1] - W_ss[0]) / (t_ss[-1] - t_ss[0])

def average_amplitude(history, method, discard_fraction = 0.75):
    q = np.array(history[method]["q"])
    N = int(len(q) * discard_fraction)
    q_ss = q[N:]
    return 0.5 * (np.max(q_ss) - np.min(q_ss))

#The initial value are in the graphics. To see hysteresis you mush change F0 to 8 and t_max = 150 and lot of point omega +100
def curves(params, omega0, method):
    """Study average dissipated power vs external frequency omega, varying beta"""
    
    betas = np.array([0.1*omega0, 0.2*omega0, 0.3*omega0, 0.5*omega0, 0.8*omega0])

    omegas = np.linspace(0.1, 3*omega0, 100)
    
    # Colormap
    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(betas)) for i in range(len(betas))]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), tight_layout = True)
    
    # 1. Main graphic: <P> vs omega for different damping parameter
    for i, beta in enumerate(betas):
        avg_powers = []
        amplitude = []
        print(f"Running => {beta/omega0}")
        for omega in omegas:
            
            osc = DrivenOscillation(q0=params.theta0, dq0=params.omega0, m=params.mass, gamma= params.mass * params.L **2 * beta, F0 = params.F0, omega=omega, t_max=params.t_max, dt=params.dt, system='nonlinear', L=params.L, F_external=params.F_external)
            model = osc.run()

            avg_powers.append(average_power(model.history, method))

            amplitude.append(average_amplitude(model.history, method))

        label = rf"$\beta = {beta/omega0:.1f}\,\omega_0$"
        ax1.plot(omegas/omega0, avg_powers, color=colors[i], lw=2.5, label=label)
        ax2.plot(omegas/omega0, amplitude, color = colors[i], lw = 2.5, label = label)

    ax1.axvline(1, color='k', ls='--', lw=1, alpha=0.7, label=r"$\omega_0$")
    ax1.set_xlabel(r"$\omega / \omega_0$")
    ax1.set_ylabel(r"$\langle P \rangle$ (W)")
    ax1.set_title("Average power vs external frequency")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.axvline(1, color='k', ls='--', lw=1, alpha=0.7, label=r"$\omega_0$")
    ax2.set_xlabel(r"$\omega / \omega_0$")
    ax2.set_ylabel("Amplitude (rad)")
    ax2.set_title("Steady-state amplitude vs external frequency")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    return betas, omegas

if __name__ == "__main__":

    # Create linear oscillator
    params = DrivenOscillationParams()
    omega_sq = np.sqrt(9.81/params.L)
    numerical_method = ['rk4', 'CN', 'Verlet']
    
    F0_list = [0, 0, 2.0]
    gamma_list = [0, 1, 2]

    #Remenber: w_d = np.sqrt(k/m - gamma/(2m)) and that must be real

    for idx, (F0_val, gamma_val) in enumerate(zip(F0_list, gamma_list)):

        print(f"Running F0={F0_val}, gamma={gamma_val}")
    
        osc = DrivenOscillation(q0=params.theta0, dq0=params.omega0, m=params.mass, gamma=gamma_val, F0=F0_val, omega=params.omega, t=params.t_max, dt=params.dt, system='nonlinear', L=params.L, F_external=params.F_external)
    
        model = osc.run() 
    
        plot_regime_summary(model.history, F0_val, gamma_val) #In this first picture we observe a nuemrical error in energy. What does happened? Sinceresly I dunno know. Because damping vibration animation works perfectly and 
        # I was studied through Euler integration and the energy more or less was conservative (whole system). However, it seems  that the numerical errors are huge due to nonlinearity of the system

    curves(params, omega_sq, method = "rk4")
    plt.show()
    
