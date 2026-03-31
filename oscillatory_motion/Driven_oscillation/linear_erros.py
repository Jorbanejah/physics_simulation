#Solo errores: RK4,CN, Verlet para distintos betas, omegas y en resonancia.
# convergencia y estabilidad para distinto dt. Mirar con la analitica

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from Driven_oscillation import DrivenOscillation

@dataclass
class DrivenOscilation_param():
    """Parameters for linear driven oscillator"""
    mass: float = 1.0
    gamma: float = 1.0
    k: float = 1.0

    omega: float = 2.0
    F0: float = 2.0
    F_external: str = 'cos'

    q0: float = 2.0
    dq0: float = 0.0

    dt: float = 0.01
    t_max: float = 50.0
    
    system: str = 'linear'

def analitical_numerical(params, numerical_method):

    linear  = DrivenOscillation(q0=params.q0, dq0=params.dq0, m=params.mass, gamma=params.gamma, F0=params.F0, omega=params.omega, t=params.t_max, dt=params.dt, system='linear', k=params.k, F_external=params.F_external)
    model = linear.run()
    history_num = model.history
    history_ana = model.analytical
    
    solution = {
        name: {
            'q' :history_num[name]['q'], 
            'v': history_num[name]['v']
            } 
        for name in numerical_method
    }

    # Add analytical and time
    solution['Analytical'] = {
        'x': history_ana['x'],
        'v': history_ana['v']
    }
    solution['t'] = history_num['rk4']['t']


    return solution

def error_ana_num(params, numerical_method):

    solution = analitical_numerical(params, numerical_method)

    error = {
        name: {"x": [], "v": []} for name in numerical_method
            }
    
    for _, num in enumerate(numerical_method):

        error[num]['x'] = np.abs(np.array(solution['Analytical']["x"]) - np.array(solution[f'{num}']['q']))
        error[num]['v'] = np.abs(np.array(solution['Analytical']["v"]) - np.array(solution[f'{num}']["v"]))
    
    return error, solution['t']


###
# -------------------- Plots ----------------
###

def plot_error_ana_num(params, numerical_method):

    cmap = plt.colormaps['viridis']
    colors = [cmap(i / len(numerical_method)) for i in range(len(numerical_method))]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)

    for ax_idx, forcing in enumerate(['cos', 'sin']):
        params.F_external = forcing
        error, t = error_ana_num(params, numerical_method)
        print(f"Computing {forcing}")
        ax = axes[ax_idx]

        for j, num in enumerate(numerical_method):
            ax.plot(t, error[num]['x'], color=colors[j], label=num)

        ax.set_title(f"External forcing: {forcing}")
        ax.set_xlabel("Time")
        ax.set_ylabel(r"$|y_{\mathrm{exact}} - y_{\mathrm{num}}|$")
        ax.set_yscale("log")
        ax.legend()

    fig.suptitle(
        rf"Numerical error for $\gamma={params.gamma}$, $\omega={params.omega}$, $k={params.k}$",
        fontsize=14
    )

def phase_portrait_error(params, numerical_method):

    cmap = plt.colormaps["viridis"]
    colors = [cmap(i / len(numerical_method)) for i in range(len(numerical_method))]

    fig, axes = plt.subplots(2, 3, figsize = (15, 7), tight_layout = True)

    for ax_idx, forcing in enumerate(['cos', 'sin']):

        params.F_external = forcing

        error, t = error_ana_num(params, numerical_method)

        for j, num in enumerate(numerical_method):
            axes[ax_idx, j].plot(error[num]["x"], error[num]["v"], color = colors[j], label = num)
            axes[ax_idx, j].set_xlabel(r"$|v_{\mathrm{exact}} - v_{\mathrm{num}}|$")
            axes[ax_idx, j].set_ylabel(r"$|y_{\mathrm{exact}} - y_{\mathrm{num}}|$")
            axes[ax_idx, j].legend()

    fig.suptitle(
        rf"Phase portrait error for $\gamma={params.gamma}$, $\omega={params.omega}$, $k={params.k}$", fontsize=14)
    
if __name__ == "__main__":

    params = DrivenOscilation_param()
    numerical_method = ['rk4', 'CN', 'Verlet']

    #plot_error_ana_num(params, numerical_method) #In this first one graphics wwe see how the verlet method doesn't fix correctly with analitcal method due to the fact that velet is for conservative system.On the other hand, the different errors between the two linear cos and sin system is so large so I want to see what happens.
    phase_portrait_error(params, numerical_method)
    plt.show()