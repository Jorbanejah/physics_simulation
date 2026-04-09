from dataclasses import dataclass
import numpy as np
from Driven_oscillation import DrivenOscillation
from manim import *

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
    t: int = 2000 #We need a huge time even 8000
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

    def extract_poincare_section(self, q_traj, dq_traj, omega, t_start=5):
        """Extract Poincaré section by sampling at driving period."""

        period = 2 * np.pi / omega
        dt = self.params.dt
        steps_per_period = int(period/dt)

        t_sample = np.arange(int(t_start/dt), len(q_traj), steps_per_period)
        
        return q_traj[t_sample], dq_traj[t_sample]

    def run_parameter_sweep(self):
        """Compute Poincaré sections for all alpha values."""
        print("Computing Poincaré sections...")
        
        for i, alpha in enumerate(self.alphas):
            progress = (i + 1) / len(self.alphas)
            bar_length = 12
            filled = int(progress * bar_length)
            bar = "█" * filled + "-" * (bar_length - filled)

            print(fr"[{bar}]  {progress*100:5.1f}%   $\alpha$ = {alpha:.2f}", end="\r", flush=True)
            
            # Update F0 = alpha
            F0 = alpha * self.params.m * self.parmas.L**2
            
            # Create oscillator with current parameters
            osc = DrivenOscillation(q0=self.params.q0, dq0=self.params.dq0, m=self.params.m, gamma=self.params.gamma, F0=F0, omega=self.params.omega, t = self.params.t, system=self.params.system, L=self.params.L, F_external = self.params.F_external)
            
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

    def store(self, filename="C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz"):
        """Store results."""
        np.savez(filename, **self.poincare_sections)
        print(f"Data stored in {filename}")

# Pre-compute data (run this first)
if __name__ == "__main__":
    system = DynamicalSystem()
    sections = system.run_parameter_sweep()
    system.store()
    data = np.load("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz", allow_pickle=True)
    rk4 = data["rk4"].item()

    for alpha, arr in rk4["q"].items():
        print(alpha, len(arr))

class PoincareSectionsScene(Scene):
    """Main animation scene showing Poincaré section evolution."""
    
    def construct(self):
        # Load pre-computed data
        data = np.load("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz", allow_pickle=True)
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

        #Create title
        title = Text("Driven pendulum equation").to_corner(UL)

        #Main equation   
        eq = MathTex(
            r"\ddot{q} + \frac{\gamma}{m L^2} \dot{q} + \frac{g}{L} \sin(q) = \frac{F_0}{m L^2} \cos(\omega t)"
        ).scale(0.9)

        self.play(FadeIn(title, shift=DOWN))

        self.play(Write(eq))

        self.wait()
        
        #Change equation
        transform_eq =  MathTex(
            r"\ddot{q} + \beta \dot{q} + \omega_0^2 \sin(q) = \alpha \cos(\omega t)"
        ).scale(0.9)
        
        self.play(
            Transform(eq, transform_eq),
        )

        self.wait()

        # Replace title
        new_title = Text("Poincaré section").to_corner(UL)
        self.play(Transform(title, new_title))
        self.wait()

        #Slide the equation away
        self.play(
            eq.animate.shift(3*RIGHT + DOWN)
        )

    def create_phase_portrait_axes(self):
        """Create phase space axes with grid and labels."""
        axes = Axes(
            x_range=[-np.pi, np.pi, np.pi/2],
            y_range=[-3, 3, 1],
            x_length=8,
            y_length=5,
            axis_config={"include_tip": False, "font_size": 24},
            tips=False
        ).shift(DOWN * 0.5)
        
        x_label = axes.get_x_axis_label(MathTex(r"q"))
        y_label = axes.get_y_axis_label(MathTex(r"\dot{q}"))
        
        grid = NumberPlane(
            x_range=[-np.pi, np.pi, np.pi/4],
            y_range=[-3, 3, 0.5],
            background_line_style={
                "stroke_opacity": 0.15, 
                "stroke_width": 1,
                "stroke_color": GREY
            }
        ).scale(0.7).shift(DOWN * 0.5)
        
        self.play(
            Create(grid),
            Create(axes),
            FadeIn(x_label),
            FadeIn(y_label),
            run_time=2
        )
        self.wait(0.5)
        
        return axes
    
    def animate_poincare_morphing(self, axes, rk4_data, alphas):
        """Animate Poincaré sections morphing through alpha values."""
        # Initial alpha (normal behavior)
        alpha0 = alphas[0]
        dots = self.make_dots(axes, rk4_data, alpha0)
        self.add(dots) #ensure dots are in the scene

        # Static label object
        alpha_label = MathTex("").scale(1.2).to_corner(UR)
        behavior_label = Text("", font_size=24).next_to(alpha_label, DOWN)

        self.add(alpha_label, behavior_label)

        # Set initial label values
        alpha_label.become(MathTex(f"\\alpha = {alpha0:.2f}").scale(1.2).to_corner(UR))
        behavior_label.become(Text("Normal", font_size=24).next_to(alpha_label, DOWN))

        self.play(
            FadeIn(dots, scale=0.8),
            FadeIn(alpha_label),
            FadeIn(behavior_label)
        )
        self.wait(1.5)

        # --- MORPH THROUGH ALPHAS ---
        for alpha in alphas[1:]:

            new_dots = self.make_dots(axes, rk4_data, alpha)

            # Behavior classification
            if alpha < 0.8:
                behavior = "Normal"
            elif alpha < 1.5:
                behavior = "Erratic"
            else:
                behavior = "Chaotic"

            # Update labels IN PLACE
            alpha_label.become(
                MathTex(f"\\alpha = {alpha:.2f}").scale(1.2).to_corner(UR)
            )
            behavior_label.become(
                Text(behavior, font_size=24).next_to(alpha_label, DOWN)
            )

            # Replace dots
            self.play(
                ReplacementTransform(dots, new_dots),
                FadeIn(alpha_label),
                FadeIn(behavior_label),
                run_time=0.6
            
            )

            dots = new_dots
            self.wait(0.3)

        self.wait(2)
    
    def make_dots(self, axes, rk4_data, alpha):
        """Create dot group for Poincaré section at given alpha."""
        q = np.array(rk4_data["q"][alpha])
        dq = np.array(rk4_data["dq"][alpha])
        
        # Limit points for performance and clarity
        q = (q + np.pi) % (2*np.pi) - np.pi

        n_points = min(800, len(q))
        indices = np.linspace(0, len(q)-1, n_points, dtype=int)

        dots = VGroup(*[
            Dot(
                axes.c2p(q[i], dq[i]), 
                radius=0.04, 
                color=BLUE
            )
            for i in indices
        ])
        
        return dots

# Usage:
# 1. Run: python poincare_section_animation.py  (computes and stores data)
# 2. Run: manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\poincare_section_animation.py PoincareSectionsScene -o my_animation.mp4
