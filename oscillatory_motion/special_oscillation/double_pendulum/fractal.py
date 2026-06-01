"""
Double Pendulum Fractal / Phase Space Diagram
==============================================

This computes the "flip time" for each initial condition (θ₁, θ₂)
with both pendulums starting from rest (ω₁ = ω₂ = 0).

Reference: FIG. 2 shows the outcome regions for the double pendulum
where angles range from -π to π.

Color scheme:
- Green: flips within 10 units
- Red: 10-100 units  
- Purple: 100-1000 units
- Blue: 1000-10000 units
- White: doesn't flip within 10000 units
- Black curve: energetically impossible (3cos(θ₁) + cos(θ₂) = 2)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Callable
import time as time_module

State = np.ndarray


@dataclass
class FractalParams:
    """Parameters for the fractal computation."""
    # Physical parameters
    g: float = 9.81
    m1: float = 1.0
    m2: float = 1.0
    L1: float = 1.0
    L2: float = 1.0
    
    # Grid parameters
    resolution: int = 500  # Grid resolution (resolution x resolution)
    theta_range: tuple[float, float] = (-np.pi, np.pi)
    
    # Time limits
    t_max_flip: float = 10000.0  # Maximum simulation time
    dt_intermediate: float = 0.01  # Output step for intermediate times
    
    # Flip thresholds (in radians)
    flip_threshold: float = np.pi
    
    # Solver tolerances
    rtol: float = 1e-8
    atol: float = 1e-10


class FlipDetector:
    """
    Event detection for pendulum flip.
    
    A flip occurs when θ crosses ±π (goes over the top).
    We use scipy's event detection for efficiency.
    """
    
    def __init__(self, threshold: float = np.pi):
        self.threshold = threshold
    
    def check_flip(self, state: State) -> bool:
        """Check if any pendulum has flipped."""
        theta1, theta2 = state[0], state[1]
        return abs(theta1) > self.threshold or abs(theta2) > self.threshold
    
    def get_flip_time(self, theta1: float, theta2: float) -> float:
        """Get the actual flip time (when crossing first occurred)."""
        # When theta = pi, that's the flip point
        return None  # Will be computed during integration


def compute_energy(state: State, params: FractalParams) -> float:
    """
    Compute total mechanical energy.
    
    Starting from rest (ω₁ = ω₂ = 0), the initial energy is all potential:
    E = -(m₁+m₂)gL₁cos(θ₁) - m₂gL₂cos(θ₂)
    
    To flip: need at least enough energy for one mass to reach the top.
    For θ₁ to flip: E ≥ (m₁+m₂)gL₁
    For θ₂ to flip: E ≥ m₂gL₂ (measured from θ₁ position)
    
    The condition for any flip is derived from energy conservation.
    """
    theta1, theta2 = state[0], state[1]
    m1, m2 = params.m1, params.m2
    L1, L2 = params.L1, params.L2
    g = params.g
    
    # Total energy (starting from rest)
    E = -(m1 + m2) * g * L1 * np.cos(theta1) - m2 * g * L2 * np.cos(theta2)
    
    return E


def can_flip(theta1: float, theta2: float, params: FractalParams) -> bool:
    """
    Check if flip is energetically possible.
    
    From energy conservation, the condition is:
    3cos(θ₁) + cos(θ₂) ≥ 2
    
    (This comes from: E ≥ (m₁+m₂)gL₁ for θ₁ flipping, with m₁=m₂=L₁=L₂)
    """
    # Critical energy for θ₁ to flip (from downward vertical)
    E_threshold = (params.m1 + params.m2) * params.g * params.L1
    
    # The blank zone equation: 3cos(θ₁) + cos(θ₂) = 2
    # Inside this curve (where 3cos(θ₁) + cos(θ₂) > 2), flipping is possible
    value = 3 * np.cos(theta1) + np.cos(theta2)
    
    return value >= 2.0


class DoublePendulumFlipSimulator:
    """Simulator optimized for flip detection."""
    
    def __init__(self, params: FractalParams):
        self.params = params
    
    def equations(self, t: float, state: State) -> State:
        """Equations of motion."""
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
        """
        Run simulation and return flip time.
        
        Returns:
            flip_time if flip occurred, -1 if no flip within t_max
        """
        params = self.params
        
        # Initial state [θ₁, θ₂, ω₁, ω₂]
        y0 = np.array([theta1_0, theta2_0, 0.0, 0.0], dtype=float)
        
        # Check if energetically possible
        if not can_flip(theta1_0, theta2_0, params):
            return -2.0  # Energetically forbidden
        
        # Time span
        t_span = (0.0, params.t_max_flip)
        
        # Use Radau for best accuracy
        sol = solve_ivp(
            fun=self.equations,
            t_span=t_span,
            y0=y0,
            method="Radau",
            max_step=params.dt_intermediate,
            rtol=params.rtol,
            atol=params.atol,
            # Event: detect when theta > pi
            events=self._create_flip_event()
        )
        
        if sol.t_events:
            # Flip detected
            return sol.t_events[0][0]
        else:
            # No flip
            return -1.0
    
    def _create_flip_event(self):
        """Create event functions for flip detection."""
        threshold = self.params.flip_threshold
        
        def event_theta1(t, y):
            return abs(y[0]) - threshold
        event_theta1.terminal = True
        event_theta1.direction = 0  # Both directions
        
        def event_theta2(t, y):
            return abs(y[1]) - threshold
        event_theta2.terminal = True
        event_theta2.direction = 0
        
        return [event_theta1, event_theta2]


def compute_fractal(params: FractalParams = None) -> tuple:
    """
    Compute the double pendulum fractal.
    
    Returns:
        theta1_grid, theta2_grid, flip_times
    """
    if params is None:
        params = FractalParams()
    
    print("=" * 60)
    print("DOUBLE PENDULUM FRACTAL COMPUTATION")
    print("=" * 60)
    print(f"Grid size: {params.resolution} x {params.resolution}")
    print(f"Time limit: {params.t_max_flip}")
    print(f"Expected simulations: {params.resolution ** 2}")
    
    start_time = time_module.time()
    
    # Create grid
    theta_min, theta_max = params.theta_range
    theta_vals = np.linspace(theta_min, theta_max, params.resolution)
    
    # Create simulator
    simulator = DoublePendulumFlipSimulator(params)
    
    # Storage
    flip_times = np.full((params.resolution, params.resolution), np.nan)
    
    # For progress tracking
    total = params.resolution ** 2
    count = 0
    last_percent = 0
    
    # Compute for each initial condition
    for i, theta1 in enumerate(theta_vals):
        for j, theta2 in enumerate(theta_vals):
            count += 1
            
            # Progress
            percent = int(100 * count / total)
            if percent != last_percent and percent % 10 == 0:
                print(f"Progress: {percent}% ({count}/{total})")
                last_percent = percent
            
            # Check energetics first
            if not can_flip(theta1, theta2, params):
                flip_times[i, j] = -2.0  # Energetically forbidden
                continue
            
            # Run simulation
            flip_time = simulator.run_simulation(theta1, theta2)
            flip_times[i, j] = flip_time
    
    elapsed = time_module.time() - start_time
    print(f"\nCompleted in {elapsed:.1f} seconds")
    
    return theta_vals, theta_vals, flip_times


def create_fractal_plot(theta1: np.ndarray, theta2: np.ndarray, 
                        flip_times: np.ndarray, save_file: str = None):
    """
    Create the fractal plot matching FIG. 2.
    """
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Define special colormap
    colors = [
        (0.0, 'white'),      # No flip (white)
        (0.0, '#00FF00'),     # Green - fast flip (<10)
        (0.3, '#00CC00'),    # Darker green
        (0.3, '#FF0000'),    # Red - medium flip (10-100)
        (0.6, '#CC0000'),    # Darker red  
        (0.6, '#800080'),     # Purple - slow flip (100-1000)
        (0.9, '#400040'),     # Darker purple
        (0.9, '#0000FF'),     # Blue - very slow (1000-10000)
        (1.0, '#0000CC'),     # Darker blue
    ]
    
    # Create custom colormap
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list('flip_cmap', colors, N=n_bins)
    
    # Log scale for visualization
    log_times = np.log10(np.clip(flip_times, 1, 1e10))
    
    # Plot
    im = ax.imshow(log_times, origin='lower', 
                    extent=[-np.pi, np.pi, -np.pi, np.pi],
                    cmap=cmap, vmin=-1, vmax=4, aspect='equal')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    
    # Custom tick labels
    tick_positions = [-1, 1, 2, 3, 4]
    tick_labels = ['<1', '10', '100', '1000', '>10000']
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    cbar.set_label('Flip Time (units of √(L/g))', fontsize=12)
    
    # Draw the energetically forbidden curve: 3cos(θ₁) + cos(θ₂) = 2
    theta1_dense = np.linspace(-np.pi, np.pi, 1000)
    
    # For the curve 3cos(θ₁) + cos(θ₂) = 2, solve for θ₂:
    # cos(θ₂) = 2 - 3cos(θ₁)
    # θ₂ = ±arccos(2 - 3cos(θ₁)) (within valid range)
    
    cos_theta2 = 2 - 3 * np.cos(theta1_dense)
    
    # Only valid where |cos(θ₂)| ≤ 1
    valid = np.abs(cos_theta2) <= 1
    theta2_curve = np.arccos(cos_theta2[valid])
    theta1_valid = theta1_dense[valid]
    
    # Plot both branches
    ax.plot(theta1_valid, theta2_curve, 'k-', linewidth=2, label='Energy limit')
    ax.plot(theta1_valid, -theta2_curve, 'k-', linewidth=2)
    
    ax.set_xlabel('θ₁ (radians)', fontsize=14)
    ax.set_ylabel('θ₂ (radians)', fontsize=14)
    ax.set_title('Double Pendulum Flip Time Fractal\n(θ₁, θ₂ initial angles, starting from rest)', 
                fontsize=14)
    
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)
    
    # Add reference lines
    ax.axhline(0, color='gray', alpha=0.3, linestyle='--')
    ax.axvline(0, color='gray', alpha=0.3, linestyle='--')
    
    if save_file:
        plt.savefig(save_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {save_file}")
    
    plt.show()
    
    return fig, ax


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == "__main__":
    # Use moderate resolution for reasonable computation time
    params = FractalParams(
        resolution=300,      # 300x300 = 90,000 simulations
        t_max_flip=10000.0,
        g=9.81,
        m1=1.0,
        m2=1.0, 
        L1=1.0,
        L2=1.0,
        rtol=1e-8,
        atol=1e-10
    )
    
    # Compute the fractal
    theta1, theta2, flip_times = compute_fractal(params)
    
    # Create the plot
    create_fractal_plot(theta1, theta2, flip_times, r'figures\\double_pendulum_fractal.png')
    
    # Print statistics
    print("\n" + "=" * 60)
    print("STATISTICS")
    print("=" * 60)
    
    n_forbidden = np.sum(flip_times < -1.5)
    n_no_flip = np.sum((flip_times < 0) & (flip_times > -1.5))
    n_fast = np.sum((flip_times > 0) & (flip_times <= 10))
    n_medium = np.sum((flip_times > 10) & (flip_times <= 100))
    n_slow = np.sum((flip_times > 100) & (flip_times <= 1000))
    n_very_slow = np.sum((flip_times > 1000) & (flip_times <= 10000))
    
    total = flip_times.size
    
    print(f"Energetically forbidden: {n_forbidden} ({100*n_forbidden/total:.1f}%)")
    print(f"No flip in 10000 units: {n_no_flip} ({100*n_no_flip/total:.1f}%)")
    print(f"Fast flip (<10): {n_fast} ({100*n_fast/total:.1f}%)")
    print(f"Medium flip (10-100): {n_medium} ({100*n_medium/total:.1f}%)")
    print(f"Slow flip (100-1000): {n_slow} ({100*n_slow/total:.1f}%)")
    print(f"Very slow flip (1000-10000): {n_very_slow} ({100*n_very_slow/total:.1f}%)")


# Run the fractal computation
params = FractalParams(resolution=200)  # Start with 200x200 for reasonable time
