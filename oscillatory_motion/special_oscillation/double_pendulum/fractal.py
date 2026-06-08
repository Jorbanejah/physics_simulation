"""
This computes the "flip time" for each initial condition (θ₁, θ₂)
with both pendulums starting from rest (ω₁ = ω₂ = 0).

Reference: FIG. 2 shows the outcome regions for the double pendulum
where angles range from -π to π.

Color scheme:

Green: flips within 10 units

Red: 10-100 units

Purple: 100-1000 units

Blue: 1000-10000 units

White: doesn't flip within 10000 units

Black curve: energetically impossible (3cos(θ₁) + cos(θ₂) = 2)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import solve_ivp
from dataclasses import dataclass
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


def compute_energy(state: State, params: FractalParams) -> float:
    """
    Compute total mechanical energy.
    
    Starting from rest (ω₁ = ω₂ = 0), the initial energy is all potential:
    E = -(m₁+m₂)gL₁cos(θ₁) - m₂gL₂cos(θ₂)
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
    """
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
        
        # Create events for flip detection
        threshold = params.flip_threshold
        
        def event_theta1_pos(t, y):
            return y[0] - threshold
        event_theta1_pos.terminal = True
        event_theta1_pos.direction = 1  # Crossing from below
        
        def event_theta1_neg(t, y):
            return y[0] + threshold
        event_theta1_neg.terminal = True
        event_theta1_neg.direction = -1  # Crossing from above
        
        def event_theta2_pos(t, y):
            return y[1] - threshold
        event_theta2_pos.terminal = True
        event_theta2_pos.direction = 1
        
        def event_theta2_neg(t, y):
            return y[1] + threshold
        event_theta2_neg.terminal = True
        event_theta2_neg.direction = -1
        
        events = [event_theta1_pos, event_theta1_neg, event_theta2_pos, event_theta2_neg]
        
        # Use Radau for best accuracy
        sol = solve_ivp(
            fun=self.equations,
            t_span=t_span,
            y0=y0,
            method="Radau",
            max_step=params.dt_intermediate,
            rtol=params.rtol,
            atol=params.atol,
            events=events
        )
        
        # Find the earliest flip time
        flip_times = []
        for t_event in sol.t_events:
            if len(t_event) > 0:
                flip_times.extend(t_event.tolist())
        
        if flip_times:
            return min(flip_times)
        else:
            # No flip
            return -1.0


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
    
    # Create a copy for visualization to avoid modifying original
    vis_data = flip_times.copy()
    
    # Handle special values: -2.0 (energetically forbidden) should show as white
    # but we need to make them distinct for the colormap
    is_forbidden = (flip_times < -1.5)
    is_no_flip = ((flip_times < 0) & (flip_times > -1.5))
    
    # For visualization: set forbidden to NaN (will show as background color)
    vis_data[is_forbidden] = np.nan
    
    # Define proper colormap with correct tuple format (position, r, g, b)
    # Position must be in ascending order
    colors = [
        (0.0, '#00FF00'),     # Green - fast flip (<10)
        (0.2, '#00CC00'),    # Darker green
        (0.2, '#FF0000'),    # Red - medium flip (10-100)
        (0.4, '#CC0000'),   # Darker red  
        (0.4, '#800080'),    # Purple - slow flip (100-1000)
        (0.6, '#400040'),    # Darker purple
        (0.6, '#0000FF'),   # Blue - very slow (1000-10000)
        (0.8, '#0000CC'),   # Darker blue
        (0.8, '#FFAA00'),   # Orange - didn't flip
        (1.0, '#FF8800'),   # Darker orange
    ]
    
    # Convert to proper format for LinearSegmentedColormap
    color_list = [c[1] for c in colors]
    position_list = [c[0] for c in colors]
    
    cmap = LinearSegmentedColormap.from_list('flip_cmap', list(zip(position_list, color_list)), N=256)
    
    # Log scale for visualization (only for positive flip times)
    vis_positive = vis_data.copy()
    vis_positive[vis_positive <= 0] = np.nan  # Treat non-flips as NaN
    
    # Log10 of flip times
    log_times = np.log10(np.clip(vis_positive, 1, 1e10))
    
    # Plot
    im = ax.imshow(log_times, origin='lower', 
                    extent=[-np.pi, np.pi, -np.pi, np.pi],
                    cmap=cmap, vmin=0, vmax=4, aspect='equal')
    
    # For energetically forbidden regions, we'll overlay them
    # Create a mask for forbidden regions
    if np.any(is_forbidden):
        forbidden_mask = np.zeros_like(flip_times, dtype=float)
        forbidden_mask[is_forbidden] = 1.0
        ax.imshow(forbidden_mask, origin='lower', 
                 extent=[-np.pi, np.pi, -np.pi, np.pi],
                 cmap='Greys', vmin=0, vmax=1, alpha=0.3, aspect='equal')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    
    # Custom tick labels
    tick_positions = [0, 1, 2, 3, 4]
    tick_labels = ['1', '10', '100', '1000', '10000']
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    cbar.set_label('Flip Time (units of √(L/g))', fontsize=12)
    
    # Draw the energetically forbidden curve: 3cos(θ₁) + cos(θ₂) = 2
    theta1_dense = np.linspace(-np.pi, np.pi, 1000)
    
    # For the curve 3cos(θ₁) + cos(θ₂) = 2, solve for θ₂:
    # cos(θ₂) = 2 - 3cos(θ₁)
    cos_theta2 = 2 - 3 * np.cos(theta1_dense)
    
    # Handle multiple branches properly
    # Valid range: -1 <= cos_theta2 <= 1
    valid = np.abs(cos_theta2) <= 1
    
    # Plot both branches using proper handling
    theta1_valid = theta1_dense[valid]
    theta2_upper = np.arccos(cos_theta2[valid])
    theta2_lower = -np.arccos(cos_theta2[valid])
    
    ax.plot(theta1_valid, theta2_upper, 'k-', linewidth=2, label='Energy limit')
    ax.plot(theta1_valid, theta2_lower, 'k-', linewidth=2)
    
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
        resolution=100,      # 100x100 = 10,000 simulations
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
    print(f"No flip in 10000 units: {n_no_flip} ({100*n_no_flip/ total :.1f}%")