import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from Driven_oscillation import DrivenOscillation
def compute_drift_wrapper(args):
    """Wrapper function for parallel processing - must be top-level"""
    params, beta, omega, omega0, dt, method = args
    return compute_drift(params, beta, omega, omega0, dt, method)

def compute_drift(params, beta, omega, omega0, dt, method='rk4'):
    """Compute energy drift for single (beta, omega, dt) combination"""
    try:
        osc = DrivenOscillation(
            q0=params.q0, dq0=params.dq0, m=params.mass,
            gamma=beta * 2 * params.mass,
            F0=params.F0, omega=omega,
            t_max=params.t_max, dt=dt,
            system='linear', k=params.k,
            F_external=params.F_external
        )
        
        model = osc.run()
        history = model.history  # Adjust this path as needed
        
        E_total = (np.array(history[method]['Ek']) + 
                  np.array(history[method]['Ep']) + 
                  np.array(history[method]['Wp_diss']) - 
                  np.array(history[method]['Wp_drive']))
        
        dE = np.max(E_total) - np.min(E_total)
        return float(np.max(np.abs(dE)))  # Return scalar
        
    except Exception as e:
        print(f"Error for beta={beta:.2f}, omega={omega:.2f}, dt={dt}: {e}")
        return np.nan

def energy_surface(params, omega0, dts=[0.01, 0.05, 0.1, 0.5], method='rk4', n_workers=None, use_parallel=True):
    """
    Efficient energy surface with dt sweep and parallelization
    """
    betas = np.linspace(omega0, 1.8*omega0, 8)
    omegas = np.linspace(-2.5*omega0, 2.5*omega0, 10)
    
    if n_workers is None:
        n_workers = min(mp.cpu_count(), 8)  # Conservative default
    
    print(f"Computing: {len(dts)} dts × {len(betas)} betas × {len(omegas)} omegas = {len(dts)*len(betas)*len(omegas)} combinations")
    print(f"Workers: {n_workers}, Parallel: {use_parallel}")
    
    drift_map = np.zeros((len(dts), len(betas), len(omegas)))
    
    # Sequential version (for debugging or single-core)
    if not use_parallel or n_workers == 1:
        print("Running sequentially...")
        for dt_idx, dt in enumerate(dts):
            for i, beta in enumerate(betas):
                for j, omega in enumerate(omegas):
                    drift_map[dt_idx, i, j] = compute_drift(params, beta, omega, omega0, dt, method)
                    if (dt_idx*len(betas)*len(omegas) + i*len(omegas) + j) % 10 == 0:
                        print(f"Progress: {100*(dt_idx*len(betas)*len(omegas) + i*len(omegas) + j)/(len(dts)*len(betas)*len(omegas)):.1f}%")
    
    # Parallel version
    # Replace the parallel execution section with:
    else:
        param_combinations = [(params, beta, omega, omega0, dt, method)
                         for dt in dts 
                         for beta in betas 
                         for omega in omegas]
    
        print("Starting parallel computation...")
        results = []
        batch_size = max(1, n_workers // 2)  # Process in smaller batches
    
        for i in range(0, len(param_combinations), batch_size):
            batch = param_combinations[i:i+batch_size]
            print(f"Processing batch {i//batch_size + 1}/{(len(param_combinations)-1)//batch_size + 1}")
        
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                batch_results = list(executor.map(compute_drift_wrapper, batch, timeout=300))
            results.extend(batch_results)
    
            # Fill results
        idx = 0
        for dt_idx in range(len(dts)):
            for i in range(len(betas)):
                for j in range(len(omegas)):
                    drift_map[dt_idx, i, j] = results[idx]
                    idx += 1
    
    # Plotting
    fig = plt.figure(figsize=(20, 12))
    
    for dt_idx, dt in enumerate(dts):
        # 3D Surface
        ax = fig.add_subplot(2, len(dts), dt_idx + 1, projection='3d')
        B, W = np.meshgrid(betas/omega0, omegas/omega0, indexing='ij')
        surf = ax.plot_surface(B, W, drift_map[dt_idx], cmap='viridis', edgecolor='none')
        ax.set_title(f'dt = {dt:.3f}')
        ax.set_xlabel(r'$\beta/\omega_0$')
        ax.set_ylabel(r'$\omega/\omega_0$')
        ax.set_zlabel(r'$\max|\Delta E|$')
        fig.colorbar(surf, ax=ax, shrink=0.6)
        
        # 2D Heatmap
        ax2 = fig.add_subplot(2, len(dts), len(dts) + dt_idx + 1)
        im = ax2.imshow(
            drift_map[dt_idx],
            extent=[omegas[0]/omega0, omegas[-1]/omega0, betas[0]/omega0, betas[-1]/omega0],
            origin='lower', aspect='auto', cmap='viridis'
        )
        ax2.set_title(f'Heatmap dt = {dt:.3f}')
        ax2.set_xlabel(r'$\omega/\omega_0$')
        ax2.set_ylabel(r'$\beta/\omega_0$')
        plt.colorbar(im, ax=ax2)
    
    plt.tight_layout()
    plt.suptitle(f'Energy Drift Analysis: {method.upper()}', y=0.95)
    plt.show()
    
    return drift_map, dts, betas/omega0, omegas/omega0

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

params = DrivenOscillationParams()
osc = DrivenOscillation(q0=params.q0, dq0=params.dq0, m=params.mass, gamma=params.gamma, F0=params.F0, omega=params.omega,t=params.t_max, dt=params.dt, system='linear',  k = params.k, F_external=params.F_external)
numerical_method = ['rk4', 'CN', 'Verlet']
omega0 = np.sqrt(params.k/params.mass)
drift_data = energy_surface(params, omega0, dts=[0.01, 0.05, 0.1, 0.5], method='rk4')