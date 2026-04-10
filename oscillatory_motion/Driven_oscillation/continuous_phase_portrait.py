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
    gamma: float = 0.1

    #External force
    F0: float = 2
    F_external: str = 'cos'
    omega: float = 2 * np.pi * 0.5

    #Time
    t: int = 800 #We need a huge time even 2000
    dt: float = 0.01

    #System
    system: str = 'nonlinear'

class DynamicalSystem:
    """
    Encapsulates numerical simulation for a sweep of the forcing amplitude.
    """

    def __init__(self):
        self.params = DrivenOscillationParameter()
        self.alphas = np.linspace(0.0, 20, 40)  # Range showing normal→chaotic
        self.methods = ["rk4"]

        # Allocate storage
        self.poincare_sections = {
            method: {"q": {}, "dq": {}, "q_full": {}, "dq_full": {}} for method in self.methods
        }

    def extract_poincare_section(self, q_traj, dq_traj, omega, t_start=400):
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

            print(rf"[{bar}]  {progress*100:5.1f}%   $\alpha$ = {alpha:.2f}", end="\r", flush=True)
            
            # Update F0 = alpha
            F0 = alpha * self.params.m * self.params.L**2
            
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
                self.poincare_sections[method]["q_full"][alpha] = q_full
                self.poincare_sections[method]["dq_full"][alpha] = dq_full
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
        print(alpha, len(arr), arr)

class PhasePortraitScene(Scene):
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
        
        # 3. Animate continuous pahse portrait

        self.show_phase_portraits(axes, rk4_data, alphas)

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

        #Slide the equation away
        self.play(
            eq.animate.to_corner(UR)
        )

        self.new_title = Text("Continuous Phase plot").to_corner(UL)
        self.play(Transform(title, self.new_title))
        
    def show_phase_portraits(self, axes, rk4_data, alphas):

        alpha_periodic = alphas[0]
        alpha_quasi    = alphas[len(alphas)//3]
        alpha_chaotic  = alphas[-1]

        curves = [
            self.make_continuous_phase_plot(axes, rk4_data, alpha_periodic),
            self.make_continuous_phase_plot(axes, rk4_data, alpha_quasi),
            self.make_continuous_phase_plot(axes, rk4_data, alpha_chaotic)
        ]

        labels = ["Periodic", "Quasi-periodic", "Chaotic"]

        for curve, label in zip(curves, labels):
            text = Text(label).to_corner(DL)
            self.play(FadeIn(text))
            self.play(Create(curve), run_time=10)
            self.wait(1.5)
            self.play(FadeOut(text), FadeOut(curve))

    def make_continuous_phase_plot(self, axes, rk4_data, alpha):

        q = np.array(rk4_data["q_full"][alpha])
        dq = np.array(rk4_data["dq_full"][alpha])

        # Wrap q
        q = (q + np.pi) % (2*np.pi) - np.pi

        # Convert to points
        pts = [axes.c2p(q[i], dq[i]) for i in range(len(q))]

        curve = VMobject()
        curve.set_points_as_corners(pts)
        curve.set_stroke(BLUE_C, 1.5, opacity=0.8)

        return curve

    def create_phase_portrait_axes(self):
        """Create global phase space axes for all alpha values."""
    
        # Global ranges based on your full sweep
        axes = Axes(
            x_range=[-np.pi, np.pi, np.pi/2],   # wrapped angle
            y_range=[-30, 30, 25],             # full dq range
            x_length=8,
            y_length=5.5,
            axis_config={"include_tip": False, "font_size": 24},
            tips=False
        ).shift(DOWN * 0.5)

        x_label = axes.get_x_axis_label(MathTex(r"\theta"))
        y_label = axes.get_y_axis_label(MathTex(r"\dot{\theta}"))

        # Grid matching the axes
        grid = NumberPlane(
            x_range=[-np.pi, np.pi, np.pi/4],
            y_range=[-30, 30, 25],
            background_line_style={
                "stroke_opacity": 0.12,
                "stroke_width": 1,
                "stroke_color": WHITE
            }
        ).match_width(axes).match_height(axes).move_to(axes)

        self.play(
            Create(grid),
            Create(axes),
            FadeIn(x_label),
            FadeIn(y_label),
            run_time=2
        )
        self.wait(0.5)

        return axes

# Usage:
# 1. Run: python poincare_section_animation.py  (computes and stores data)
# 2. Run: manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\continuous_phase_portrait.py PhasePortraitScene -o my_animation.mp4