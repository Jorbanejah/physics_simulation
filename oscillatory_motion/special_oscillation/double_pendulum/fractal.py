"""
Double Pendulum Fractal / Phase Space Diagram
==============================================

This computes the "flip time" for each initial condition (θ₁, θ₂)
with both pendulums starting from rest (ω₁ = ω₂ = 0).

Reference: shows the outcome regions for the double pendulum
where angles range from -π to π.

Color scheme:

    Extended colormap with many more color categories.
    
    Categories:
    Continuous viridis colormap
    Black curve: energetically impossible (3cos(θ₁) + cos(θ₂) = 2)

==============================================
Note: friendly reminder

In FractalParams - these parameters control computation time

resolution: int = 500      # ↑ Increases quadratically (N² simulations)
t_max_flip: float = 10000.0  # ↑ More time = more integration steps
dt_intermediate: float = 0.01  # ↓ Smaller step = more accurate but slower
rtol: float = 1e-8         # ↓ Tighter tolerance = slower but more accurate
atol: float = 1e-10

The current code takes about 27 hours and 43 minutes to run. So... think twice ;)


Jorge Orbaneja Huerta
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.integrate import solve_ivp
from dataclasses import dataclass
import time as time_module

State = np.ndarray


# ============================================================
# PARAMETERS
# ============================================================

@dataclass
class FractalParams:
    g: float = 9.81
    m1: float = 1.0
    m2: float = 1.0
    L1: float = 1.0
    L2: float = 1.0
    resolution: int = 200
    theta_range: tuple[float, float] = (-np.pi, np.pi)
    t_max_flip: float = 500.0
    dt_intermediate: float = 0.05
    flip_threshold: float = np.pi
    rtol: float = 1e-4
    atol: float = 1e-6


# ============================================================
# DOUBLE PENDULUM DYNAMICS
# ============================================================

class DoublePendulumFlipSimulator:
    def __init__(self, params: FractalParams):
        self.params = params
    
    def equations(self, t: float, state: State) -> State:
        theta1, theta2, omega1, omega2 = state
        m1, m2 = self.params.m1, self.params.m2
        L1, L2 = self.params.L1, self.params.L2
        g = self.params.g
        
        delta = theta1 - theta2
        sin_delta = np.sin(delta)
        cos_delta = np.cos(delta)
        
        M11 = (m1 + m2) * L1**2
        M12 = m2 * L1 * L2 * cos_delta
        M22 = m2 * L2**2
        
        F1 = -(m1 + m2) * g * L1 * np.sin(theta1) - m2 * L1 * L2 * omega2**2 * sin_delta
        F2 = m2 * L1 * L2 * omega1**2 * sin_delta - m2 * g * L2 * np.sin(theta2)
        
        det = M11 * M22 - M12**2
        alpha1 = (M22 * F1 - M12 * F2) / det
        alpha2 = (-M12 * F1 + M11 * F2) / det
        
        return np.array([omega1, omega2, alpha1, alpha2])
    
    def run_simulation(self, theta1_0: float, theta2_0: float) -> float:
        params = self.params
        y0 = np.array([theta1_0, theta2_0, 0.0, 0.0], dtype=float)

        sol = solve_ivp(
            fun=self.equations,
            t_span=(0.0, params.t_max_flip),
            y0=y0,
            method="DOP",
            max_step=params.dt_intermediate,
            rtol=params.rtol,
            atol=params.atol,
            events=self._create_flip_event()
        )
        
        # Default: no flip detected
        flip_time = np.nan
        
        if sol.t_events is not None:
            for event_times in sol.t_events:
                if len(event_times) > 0:
                    flip_time = event_times[0]
                    break
        
        return flip_time
    
    def _create_flip_event(self):
        threshold = self.params.flip_threshold
        
        def event_theta1(t, y):
            return abs(y[0]) - threshold
        event_theta1.terminal = True
        event_theta1.direction = 0
        
        def event_theta2(t, y):
            return abs(y[1]) - threshold
        event_theta2.terminal = True
        event_theta2.direction = 0
        
        return [event_theta1, event_theta2]


# ============================================================
# FRACTAL COMPUTATION
# ============================================================

def compute_fractal(params: FractalParams = None) -> tuple:
    if params is None:
        params = FractalParams()
    
    print("=" * 60)
    print("DOUBLE PENDULUM FRACTAL COMPUTATION")
    print("=" * 60)
    print(f"Grid size: {params.resolution} x {params.resolution}")
    print(f"Time limit: {params.t_max_flip}")
    
    start_time = time_module.time()
    
    theta_vals = np.linspace(*params.theta_range, params.resolution)
    simulator = DoublePendulumFlipSimulator(params)
    
    flip_times = np.full((params.resolution, params.resolution), np.nan)
    
    total = params.resolution ** 2
    count = 0
    last_percent = 0
    
    for i, theta1 in enumerate(theta_vals):
        for j, theta2 in enumerate(theta_vals):
            count += 1
            
            percent = int(100 * count / total)
            if percent != last_percent and percent % 10 == 0:
                elapsed = time_module.time() - start_time
                eta = elapsed * (100 - percent) / percent if percent > 0 else 0
                print(f"Progress: {percent}% ({count}/{total}) - ETA: {eta:.0f}s")
                last_percent = percent
            
            flip_time = simulator.run_simulation(theta1, theta2)
            flip_times[j, i] = flip_time
    
    elapsed = time_module.time() - start_time
    print(f"\nCompleted in {elapsed:.1f} seconds")
    
    return theta_vals, theta_vals, flip_times


# ============================================================
# FRACTAL PLOT (VIRIDIS)
# ============================================================

def create_fractal_plot(theta1, theta2, flip_times, save_file=None):
    fig, ax = plt.subplots(figsize=(12, 10))

    # Continuous viridis colormap
    display_data = flip_times.copy()
    display_data[display_data <= 0] = np.nan

    im = ax.imshow(
        display_data,
        origin='lower',
        extent=[-np.pi, np.pi, -np.pi, np.pi],
        cmap='viridis',
        aspect='equal'
    )

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Flip Time', fontsize=12)

    ax.set_xlabel('θ₁ (radians)', fontsize=14)
    ax.set_ylabel('θ₂ (radians)', fontsize=14)
    ax.set_title('Double Pendulum Flip Time Fractal', fontsize=14)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)

    if save_file:
        import os
        directory = os.getcwd()
        route = os.path.join(directory, save_file)
        plt.savefig(route, dpi=150, bbox_inches='tight')
        print(f"Saved to {route}")

    plt.show()
    return fig, ax


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    params = FractalParams(
        resolution=200,
        t_max_flip=300.0,
        rtol=1e-4,
        atol=1e-6
    )
    
    theta1, theta2, flip_times = compute_fractal(params)
    create_fractal_plot(theta1, theta2, flip_times, save_file='fractal_viridis.png')