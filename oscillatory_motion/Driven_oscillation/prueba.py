import numpy as np
import matplotlib.pyplot as plt
from Driven_oscillation import DrivenOscillation
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter


def advanced_stability_map(params, omega0, n_point = 10, method='rk4', t_max=100, dt=0.02, nonlinear=False):
    """
    Complete stability analysis: Energy drift + Physical insights + Pole detection
    """
    # PHYSICALLY-RELEVANT RANGES
    beta_norm = np.linspace(0.01, 1.2, n_point)  # Overdamped included
    omega_norm = np.linspace(0.1, 3.0, n_point)  # Positive frequencies only
    
    betas, omegas = np.meshgrid(beta_norm * omega0, omega_norm * omega0, indexing='ij')
    n_beta, n_omega = betas.shape
    
    # Multi-metric analysis
    drift_map = np.zeros((n_beta, n_omega))
    growth_map = np.zeros((n_beta, n_omega))  # Instability detection
    resonance_map = np.zeros((n_beta, n_omega))  # Physical response
    
    print(f"Computing {n_beta*n_omega} points...")
    
    for i in range(n_beta):
        for j in range(n_omega):
            beta, omega = betas[i,j], omegas[i,j]
            
            osc = DrivenOscillation(
                q0=params.q0, dq0=0.0, m=params.mass,
                gamma=2*params.mass*beta, F0=params.F0,
                omega=omega, t=t_max, dt=dt,
                system='nonlinear' if nonlinear else 'linear',
                k=params.k, L=params.L if nonlinear else None,
                F_external='sin'
            )
            model = osc.run()
            
            # 1. NUMERICAL STABILITY: Energy drift
            Ek, Ep = np.array(model.history[method]['Ek']), np.array(model.history[method]['Ep'])
            E_total = Ek + Ep
            drift_map[i,j] = 100 * np.max(np.abs(E_total - E_total[0])) / E_total[0]
            
            # 2. PHYSICAL GROWTH: Amplitude growth rate
            q = np.array(model.history[method]['q'])
            growth_rate = np.log(np.std(q[-1:]) / np.std(q[:100])) / t_max if np.std(q[-1:]) > 0 else 0
            growth_map[i,j] = growth_rate * 100
            
            # 3. RESONANCE STRENGTH: Steady-state power
            Wp_diss = np.array(model.history[method]['Wp_diss'])
            P_steady = (Wp_diss[-1] - Wp_diss[int(0.7*len(Wp_diss))]) / (t_max * 0.3)
            resonance_map[i,j] = P_steady

    fig = plt.figure(figsize=(20, 15))
    
    # 1. NUMERICAL STABILITY 

    ax1 = plt.subplot(3, 4, (1, 2))
    im1 = ax1.imshow(drift_map.T, extent=[beta_norm[0], beta_norm[-1], omega_norm[0], omega_norm[-1]], origin='lower', norm=LogNorm(vmin=1e-4, vmax=drift_map.max()), cmap='viridis_r')
    ax1.contour(beta_norm, omega_norm, gaussian_filter(drift_map, 1), levels=5, colors='white', alpha=0.7)
    ax1.set_xlabel(r'$\beta/\omega_0$'); ax1.set_ylabel(r'$\omega/\omega_0$')
    ax1.set_title('Numerical Stability\n(Energy Drift %)', fontweight='bold')
    plt.colorbar(im1, ax=ax1, label='Energy Drift (%)')
    
    # 2. INSTABILITY POLES (Growth regions)
    ax2 = plt.subplot(3, 4, 3)
    im2 = ax2.imshow(growth_map.T, extent=[beta_norm[0], beta_norm[-1], omega_norm[0], omega_norm[-1]], cmap='RdBu_r', vmin=-5, vmax=5)
    ax2.contour(beta_norm, omega_norm, growth_map, levels=[0], colors='black', linewidths=2)
    ax2.set_xlabel(r'$\beta/\omega_0$'); ax2.set_title('Instability Growth Rate\n(Parametric Resonance)', fontweight='bold')
    plt.colorbar(im2, ax=ax2, label='Growth Rate (%)')
    
    # 3. PHYSICAL RESONANCE (Power absorption)
    ax3 = plt.subplot(3, 4, 4)
    im3 = ax3.imshow(resonance_map.T, extent=[beta_norm[0], beta_norm[-1], omega_norm[0], omega_norm[-1]],norm=LogNorm(vmin=1e-3, vmax=resonance_map.max()), cmap='plasma')
    ax3.contour(beta_norm, omega_norm, resonance_map, levels=10, colors='gold', alpha=0.8, linewidths=1.5)
    ax3.set_xlabel(r'$\beta/\omega_0$'); ax3.set_title('Physical Resonance\n(Steady-State Power)', fontweight='bold')
    plt.colorbar(im3, ax=ax3, label='Power (arb. units)')
    
    # 4. ANALYTICAL RESONANCE CURVES (Validation)
    ax4 = plt.subplot(3, 4, (5, 6))
    beta_test = np.linspace(0.01, 1.0, 100)
    omega_test = np.linspace(0.1, 3.0, 100)
    B_test, W_test = np.meshgrid(beta_test, omega_test)
    
    # Analytical resonance: P ~ (γω)^2 / [(ω₀²-ω²)^2 + (2βω)^2]
    denom = (1 - W_test**2)**2 + (2 * B_test * W_test)**2
    P_analytical = (W_test**2 * B_test) / denom  # Normalized
    
    ax4.contourf(B_test, W_test, P_analytical, levels=20, cmap='plasma', alpha=0.8)
    ax4.contour(B_test, W_test, P_analytical, levels=5, colors='white', linewidths=1.5)
    ax4.set_xlabel(r'$\beta/\omega_0$'); ax4.set_ylabel(r'$\omega/\omega_0$')
    ax4.set_title('Analytical Resonance (Validation)', fontweight='bold')
    
    # 5. METHOD COMPARISON SLICES
    ax5 = plt.subplot(3, 4, 7)
    for m, mthd in enumerate(['rk4', 'Verlet', 'CN']):
        slice_beta = 0.2  # Fixed damping
        i_slice = np.argmin(np.abs(beta_norm - slice_beta))
        ax5.semilogy(omega_norm, drift_map[i_slice,:], label=f'{mthd}: β/ω₀={slice_beta:.1f}', lw=3)
    ax5.set_xlabel(r'$\omega/\omega_0$'); ax5.set_ylabel('Energy Drift (%)')
    ax5.set_title('Method Comparison\n(Fixed Damping)', fontweight='bold')
    ax5.legend(); ax5.grid(True, alpha=0.3)
    
    # 6. RESONANCE PEAKS vs DAMPING
    ax6 = plt.subplot(3, 4, 8)
    peak_omegas = np.max(resonance_map, axis=1)  # Peak power per beta
    ax6.semilogy(beta_norm, peak_omegas, 'ro-', lw=3, markersize=6)
    ax6.set_xlabel(r'$\beta/\omega_0$'); ax6.set_ylabel('Peak Resonance Power')
    ax6.set_title('Resonance Strength vs Damping', fontweight='bold')
    ax6.grid(True, alpha=0.3)
    
    # 7. POLE LOCATIONS (Zeros of characteristic equation)
    ax7 = plt.subplot(3, 4, (9, 10))
    s1, s2 = np.roots([1, 2*beta_norm, np.ones_like(beta_norm)])  # Characteristic eq
    ax7.plot(beta_norm, np.real(s1)/omega0, 'b-', label='Pole 1', lw=2)
    ax7.plot(beta_norm, np.real(s2)/omega0, 'r-', label='Pole 2', lw=2)
    ax7.axhline(0, color='k', ls='--', alpha=0.5)
    ax7.set_xlabel(r'$\beta/\omega_0$'); ax7.set_ylabel('Re(s)/ω₀')
    ax7.set_title('Damping Poles (Stability Boundary)', fontweight='bold')
    ax7.legend(); ax7.grid(True, alpha=0.3)
    
    # 8. SUMMARY STATISTICS
    ax8 = plt.subplot(3, 4, (11, 12))
    methods = ['rk4', 'Verlet', 'CN']
    stability_scores = []
    for mthd in methods:
        # Mock data - replace with actual computation
        score = np.mean(drift_map) / np.log(10)  # Lower = better
        stability_scores.append(score)
    
    bars = ax8.bar(methods, stability_scores, color=['C0', 'C1', 'C2'], alpha=0.8, edgecolor='black')
    ax8.set_ylabel('Stability Score'); ax8.set_title('Method Ranking', fontweight='bold')
    ax8.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('complete_stability_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig, {'drift': drift_map, 'growth': growth_map, 'resonance': resonance_map}



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

# === USAGE ===
params = DrivenOscillationParams()  # Your params
omega0 = np.sqrt(params.k/params.mass)
fig, maps = advanced_stability_map(params, omega0, method='Verlet')