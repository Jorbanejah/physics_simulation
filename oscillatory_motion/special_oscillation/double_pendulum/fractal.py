"""
Double Pendulum Fractal / Phase Space Diagram
==============================================

This computes the "flip time" for each initial condition (θ₁, θ₂)
with both pendulums starting from rest (ω₁ = ω₂ = 0).

Reference: FIG. 2 shows the outcome regions for the double pendulum
where angles range from -π to π.

Color scheme:

    Extended colormap with many more color categories.
    
    Categories:
    0: Forbidden (white)
    1: No flip (light gray)
    2: Very fast < 1 (bright green)
    3: Fast 1-3 (green)  
    4: Medium-fast 3-10 (lime)
    5: Medium 10-30 (yellow)
    6: Medium-slow 30-100 (yellow-orange)
    7: Slow 100-300 (orange)
    8: Very slow 300-1000 (red-orange)
    9: Extremely slow > 1000 (red)
    10: Black curve: energetically impossible (3cos(θ₁) + cos(θ₂) = 2)

Double Pendulum Fractal / Phase Space Diagram
==============================================
Optimized version with faster solver and early termination
"""
"""
Extended colormap version with correct color bins
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from scipy.integrate import solve_ivp
from dataclasses import dataclass
import time as time_module

State = np.ndarray


@dataclass
class FractalParams:
    """Parameters for the fractal computation."""
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


def can_flip(theta1: float, theta2: float, params: FractalParams) -> bool:
    value = 3 * np.cos(theta1) + np.cos(theta2)
    return value <=-2


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
        
        if not can_flip(theta1_0, theta2_0, params):
            return -2.0
        
        t_span = (0.0, params.t_max_flip)
        
        sol = solve_ivp(
            fun=self.equations,
            t_span=t_span,
            y0=y0,
            method="RK45",
            max_step=params.dt_intermediate,
            rtol=params.rtol,
            atol=params.atol,
            events=self._create_flip_event()
        )
        
        flip_time = -1.0
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


def compute_fractal(params: FractalParams = None) -> tuple:
    if params is None:
        params = FractalParams()
    
    print("=" * 60)
    print("DOUBLE PENDULUM FRACTAL COMPUTATION")
    print("=" * 60)
    print(f"Grid size: {params.resolution} x {params.resolution}")
    print(f"Time limit: {params.t_max_flip}")
    print(f"Expected simulations: {params.resolution ** 2}")
    
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
            
            if not can_flip(theta1, theta2, params):
                flip_times[i, j] = -2.0
                continue
            
            flip_time = simulator.run_simulation(theta1, theta2)
            flip_times[j, i] = flip_time
    
    elapsed = time_module.time() - start_time
    print(f"\nCompleted in {elapsed:.1f} seconds")
    
    return theta_vals, theta_vals, flip_times


def create_fractal_plot_extended(theta1, theta2, flip_times, 
                                 use_discrete=True, 
                                 save_file=None):
    """Create the fractal plot with extended colormap."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    if use_discrete:
        # Extended discrete colormap - CORRECTED bin counts
        categories = np.full_like(flip_times, np.nan)
        categories[flip_times < -1.5] = 0          # Forbidden
        categories[(flip_times >= -1.5) & (flip_times < 0)] = 1  # No flip
        categories[(flip_times > 0) & (flip_times <= 1)] = 2     # Very fast
        categories[(flip_times > 1) & (flip_times <= 3)] = 3    # Fast
        categories[(flip_times > 3) & (flip_times <= 10)] = 4     # Medium-fast
        categories[(flip_times > 10) & (flip_times <= 30)] = 5     # Medium
        categories[(flip_times > 30) & (flip_times <= 100)] = 6  # Medium-slow
        categories[(flip_times > 100) & (flip_times <= 300)] = 7  # Slow
        categories[(flip_times > 300) & (flip_times <= 1000)] = 8   # Very slow
        categories[flip_times > 1000] = 9                         # Extremely slow
        
        colors = [
            '#FFFFFF',  # 0: Forbidden (white)
            '#D0D0D0',  # 1: No flip (gray)
            '#00FF00',  # 2: Very fast (bright green)
            '#00DD00',  # 3: Fast (green)
            '#44FF44',  # 4: Medium-fast (lime)
            '#CCFF00',  # 5: Medium (yellow-green)
            '#FFCC00',  # 6: Medium-slow (gold)
            '#FF8800',  # 7: Slow (orange)
            '#FF4400',  # 8: Very slow (red-orange)
            '#FF0000',  # 9: Extremely slow (red)
        ]
        
        cmap = ListedColormap(colors)
        bounds = np.arange(len(colors) + 1)  # [0, 1, 2, ..., 10]
        norm = BoundaryNorm(bounds, cmap.N)
        
        im = ax.imshow(categories, origin='lower', 
                       extent=[-np.pi, np.pi, -np.pi, np.pi],
                       cmap=cmap, norm=norm, aspect='equal')
        
        # Custom colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8, 
                          ticks=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
        cbar.set_ticklabels([
            'Forbidden',
            'No flip',
            '0-1',
            '1-3',
            '3-10',
            '10-30',
            '30-100',
            '100-300',
            '300-1000',
            '>1000'
        ])
        cbar.set_label('Flip Time', fontsize=12)
        
    else:
        # Smooth continuous colormap
        display_data = flip_times.copy()
        display_data[display_data <= 0] = np.nan
        log_data = np.log10(display_data + 0.1)
        
        # Smooth gradient colors
        colors_list = [
            (0.0, '#FFFFFF'),   # White (forbidden/no flip)
            (0.01, '#00FF00'),  # Green
            (0.15, '#44FF44'),  # Lime
            (0.25, '#AAFF00'),  # Yellow-green
            (0.35, '#FFFF00'),  # Yellow
            (0.45, '#FFAA00'),  # Gold
            (0.55, '#FF6600'),  # Orange
            (0.65, '#FF3300'),  # Red-orange
            (0.75, '#FF0000'),  # Red
            (0.85, '#AA00FF'),  # Purple
            (0.95, '#4400FF'),  # Blue
            (1.0, '#0000AA'),   # Dark blue
        ]
        
        cmap = LinearSegmentedColormap.from_list('fractal_cmap', colors_list, N=256)
        
        im = ax.imshow(log_data, origin='lower', 
                       extent=[-np.pi, np.pi, -np.pi, np.pi],
                       cmap=cmap, vmin=-1, vmax=3, aspect='equal')
        
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        tick_vals = [-1, 0, 1, 2, 3]
        tick_labels = ['0.1', '1', '10', '100', '1000']
        cbar.set_ticks(tick_vals)
        cbar.set_ticklabels(tick_labels)
        cbar.set_label('Flip Time (log scale)', fontsize=12)
    
    # Draw energetically forbidden curve
    theta1_dense = np.linspace(-np.pi, np.pi, 1000)
    cos_theta2 = 2 - 3 * np.cos(theta1_dense)
    valid = np.abs(cos_theta2) <= 1
    theta2_curve = np.arccos(cos_theta2[valid])
    theta1_valid = theta1_dense[valid]
    
    ax.plot(theta1_valid, theta2_curve, 'k-', linewidth=2)
    ax.plot(theta1_valid, -theta2_curve, 'k-', linewidth=2)
    
    ax.set_xlabel('θ₁ (radians)', fontsize=14)
    ax.set_ylabel('θ₂ (radians)', fontsize=14)
    ax.set_title('Double Pendulum Flip Time Fractal\n(Extended Colormap)', fontsize=14)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)
    
    ax.axhline(0, color='gray', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.axvline(0, color='gray', alpha=0.3, linestyle='--', linewidth=0.5)
    
    if save_file:
        plt.savefig(save_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {save_file}")
    
    plt.show()
    return fig, ax


if __name__ == "__main__":
    params = FractalParams(
        resolution=10,
        t_max_flip=100000.0,
        rtol=1e-4,
        atol=1e-6
    )
    
    theta1, theta2, flip_times = compute_fractal(params)
    create_fractal_plot_extended(theta1, theta2, flip_times, 
                                use_discrete=True,
                                save_file='fractal_extended.png')
"""
Note: friendly reminder

In FractalParams - these parameters control computation time

resolution: int = 500      # ↑ Increases quadratically (N² simulations)
t_max_flip: float = 10000.0  # ↑ More time = more integration steps
dt_intermediate: float = 0.01  # ↓ Smaller step = more accurate but slower
rtol: float = 1e-8         # ↓ Tighter tolerance = slower but more accurate
atol: float = 1e-10

"""