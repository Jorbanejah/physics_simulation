from manim import *
import numpy as np

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