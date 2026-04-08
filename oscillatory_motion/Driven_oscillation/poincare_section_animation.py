from dataclasses import dataclass
import numpy as np
from Driven_oscillation import DrivenOscillation
import manim as m

@dataclass
class DrivenOscillationParameter:
    #Initial condition
    q0: float = np.deg2rad(30)
    dq0: float = 0

    #Innate parameter
    L: float = 1
    m: float = 2
    gamma: float = 0.4

    #External force
    F0: float = 2
    F_external: str = 'cos'
    omega: float = 2 * np.pi * 0.8

    #Time
    t: int = 2000
    dt: float = 0.01

    #System
    system: str = 'nonlinear'

class DynamicalSystem:
    """
    Encapsulates numerical simulation for a sweep of the forcing amplitude.
    """

    def __init__(self):
        self.params = DrivenOscillationParameter()
        self.alphas = np.linspace(0.3, 2.5, 30)  # Range showing normal→chaotic
        self.methods = ["rk4", "CN", "Verlet"]

        # Allocate storage
        self.poincare_sections = {
            method: {"q": {}, "dq": {}} for method in self.methods
        }

    def extract_poincare_section(self, q_traj, dq_traj, omega, t_start=500):
        """Extract Poincaré section by sampling at driving period."""

        period = 2 * np.pi / omega
        dt = self.params.dt
        t_sample = np.arange(t_start, len(q_traj) * dt, period)
        t_sample = t_sample[t_sample/dt < len(q_traj)].astype(int)
        
        return q_traj[t_sample], dq_traj[t_sample]

    def run_parameter_sweep(self):
        """Compute Poincaré sections for all alpha values."""
        print("Computing Poincaré sections...")
        
        for i, alpha in enumerate(self.alphas):
            print(f"Alpha {i+1}/{len(self.alphas)}: {alpha:.2f}")
            
            # Update F0 = alpha
            F0 = alpha
            
            # Create oscillator with current parameters
            osc = DrivenOscillation(q0=self.params.q0, dq0=self.params.dq0, m=self.params.m, gamma=self.params.gamma, F0=F0, omega=self.params.omega, system=self.params.system, L=self.params.L, F_external = self.params.F_external)
            
            # Run simulation
            model = osc.run()
            
            # Extract Poincaré sections for each method
            for method in self.methods:
                q_full = np.array(model.history[method]["q"])
                dq_full = np.array(model.history[method]["v"])
                
                # Extract Poincaré section (discard transients)
                q_poincare, dq_poincare = self.extract_poincare_section(q_full, dq_full, self.params.omega)
                
                self.poincare_sections[method]["q"][alpha] = q_poincare
                self.poincare_sections[method]["dq"][alpha] = dq_poincare
        
        return self.poincare_sections

    def store(self, filename="poincare_sections.npz"):
        """Store results."""
        np.savez(filename, **self.poincare_sections)
        print(f"Data stored in {filename}")

# Pre-compute data (run this first)
if __name__ == "__main__":
    system = DynamicalSystem()
    sections = system.run_parameter_sweep()
    system.store()

class PoincareSectionsScene(m.Scene):
    """Main animation scene showing Poincaré section evolution."""
    
    def construct(self):
        # Load pre-computed data
        data = np.load("poincare_sections.npz", allow_pickle=True)
        rk4_data = data["rk4"].item()
        
        alphas = sorted(rk4_data["q"].keys())
        
        # 1. Show equation
        self.show_equation()
        
        # 2. Create axes, labels, and grid
        axes = self.create_phase_portrait_axes()
        
        # 3. Animate Poincaré sections morphing through alpha values
        self.animate_poincare_morphing(axes, rk4_data, alphas)
    
    def show_equation(self):
        """Display the driven pendulum equation."""
        eq = m.MathTex(
            r"\ddot{q} + \gamma \dot{q} + \frac{g}{L}\sin(q) = \frac{F_0}{m}\cos(\omega t)"
        ).scale(0.9).to_edge(m.UP)
        
        self.play(m.Write(eq))
        self.wait(1.5)
    
    def create_phase_portrait_axes(self):
        """Create phase space axes with grid and labels."""
        axes = m.Axes(
            x_range=[-np.pi, np.pi, np.pi/2],
            y_range=[-3, 3, 1],
            x_length=8,
            y_length=5,
            axis_config={"include_tip": False, "font_size": 24},
            tips=False
        ).shift(m.DOWN * 0.5)
        
        x_label = axes.get_x_axis_label(m.MathTex(r"q"))
        y_label = axes.get_y_axis_label(m.MathTex(r"\dot{q}"))
        
        grid = m.NumberPlane(
            x_range=[-np.pi, np.pi, np.pi/4],
            y_range=[-3, 3, 0.5],
            background_line_style={
                "stroke_opacity": 0.15, 
                "stroke_width": 1,
                "stroke_color": m.GREY
            }
        ).scale(0.7).shift(m.DOWN * 0.5)
        
        self.play(
            m.Create(grid),
            m.Create(axes),
            m.FadeIn(x_label),
            m.FadeIn(y_label),
            run_time=2
        )
        self.wait(0.5)
        
        return axes
    
    def animate_poincare_morphing(self, axes, rk4_data, alphas):
        """Animate Poincaré sections morphing through alpha values."""
        # Initial alpha (normal behavior)
        alpha0 = alphas[0]
        dots = self.make_dots(axes, rk4_data, alpha0)
        
        alpha_label = m.MathTex(f"\\alpha = {alpha0:.2f}").scale(1.2).to_corner(m.UR)
        behavior_label = m.Text("Normal", font_size=24).next_to(alpha_label, m.DOWN)
        
        self.play(
            m.FadeIn(dots, scale=0.8),
            m.Write(alpha_label),
            m.Write(behavior_label)
        )
        self.wait(1.5)
        
        # Morph through alphas
        for i, alpha in enumerate(alphas[1:]):
            new_dots = self.make_dots(axes, rk4_data, alpha)
            
            # Update behavior label based on alpha
            if alpha < 0.8:
                behavior = "Normal"
                color = m.BLUE
            elif alpha < 1.5:
                behavior = "Erratic"
                color = m.YELLOW
            else:
                behavior = "Chaotic"
                color = m.RED
            
            new_alpha_label = m.MathTex(f"\\alpha = {alpha:.2f}").scale(1.2).to_corner(m.UR)
            new_behavior_label = m.Text(behavior, font_size=24).next_to(new_alpha_label, m.DOWN)
            
            self.play(
                m.Transform(dots, new_dots),
                m.Transform(alpha_label, new_alpha_label),
                m.Transform(behavior_label, new_behavior_label),
                run_time=0.6
            )
            
            dots = new_dots
            alpha_label = new_alpha_label
            behavior_label = new_behavior_label
            
            self.wait(0.3)
        
        self.wait(2)
    
    def make_dots(self, axes, rk4_data, alpha):
        """Create dot group for Poincaré section at given alpha."""
        q = np.array(rk4_data["q"][alpha])
        dq = np.array(rk4_data["dq"][alpha])
        
        # Limit points for performance and clarity
        n_points = min(800, len(q))
        indices = np.linspace(0, len(q)-1, n_points, dtype=int)
        
        dots = m.VGroup(*[
            m.Dot(
                axes.c2p(q[idx], dq[idx]), 
                radius=0.035, 
                color=m.BLUE
            )
            for idx in indices
        ])
        
        return dots

# Usage:
# 1. Run: python poincare_section_animation.py  (computes and stores data)
# 2. Run: manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\poincare_section_animation.py PoincareSectionsScene -o my_animation.mp4
