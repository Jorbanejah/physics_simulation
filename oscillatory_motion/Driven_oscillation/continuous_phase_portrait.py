from manim import *
import numpy as np

class PhasePortraitScene(Scene):
    """Main animation scene showing Poincaré section evolution."""
    
    def construct(self):
        # Load pre-computed data (adjust path as needed)
        data = np.load("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz", allow_pickle=True)
        
        q_poincare = data["poincare_q"].item()
        dq_poincare = data["poincare_dq"].item()

        y = q_poincare, dq_poincare
        alphas = sorted(dq_poincare.keys())
        
        # 1. Show equation
        self.show_equation()
        
        # 2. Create axes, labels, and grid
        axes = self.create_phase_portrait_axes()
        
        # 3. Animate Poincaré sections progressively
        self.show_poincare_sections(axes, y, alphas)

    def show_equation(self):
        """Display the driven pendulum equation."""
        # Create title
        title = Text("Driven Pendulum", font_size=36).to_corner(UL)

        # Main equation   
        eq = MathTex(
            r"\ddot{\theta} + \gamma\,\dot{\theta} + \sin(\theta) = f\,\cos(\omega t)"
        ).scale(0.9).next_to(title, DOWN, buff=0.5)

        self.play(FadeIn(title), Write(eq), run_time=2)
        self.wait(1)
        
        # Slide equation away
        self.play(eq.animate.shift(RIGHT*6 + DOWN*0.5))
        poincare_title = Text("Poincaré Sections", font_size=36).to_corner(UL)
        self.play(Transform(title, poincare_title))
        
    def create_phase_portrait_axes(self):
        """Create phase space axes for Poincaré sections."""
        axes = Axes(
            x_range=[-np.pi, np.pi, np.pi/2],
            y_range=[-3.5, 3.5, 1],
            x_length=8,
            y_length=4.5,
            axis_config={"include_tip": False, "font_size": 24},
            tips=False
        ).shift(DOWN * 0.3)

        x_label = axes.get_x_axis_label(MathTex(r"\theta"))
        y_label = axes.get_y_axis_label(MathTex(r"\dot{\theta}"))

        grid = NumberPlane(
            x_range=axes.x_range,
            y_range=axes.y_range,
            x_length=8,
            y_length=4.5,
            background_line_style={
                "stroke_opacity": 0.15,
                "stroke_width": 1,
                "stroke_color": WHITE
            }
        ).move_to(axes)

        section_line = Line(
            axes.c2p(-np.pi, 0), axes.c2p(np.pi, 0),
            color=YELLOW, stroke_width=4, stroke_opacity=0.8
        )
        section_label = Text("Poincaré Section", font_size=24, color=YELLOW).next_to(section_line, DOWN, buff=0.2)

        self.play(
            Create(grid),
            Create(axes),
            FadeIn(x_label),
            FadeIn(y_label),
            Create(section_line),
            Write(section_label),
            run_time=2
        )
        self.wait(0.5)

        return axes
    
    def show_poincare_sections(self, axes, y, alphas):
        """Animate Poincaré sections building up progressively."""
        theta_dict, omega_dict = y
        
        # Color cycle for different regimes
        colors = [BLUE_C, GREEN_C, ORANGE, RED_C]
        
        all_points = VGroup()
        
        for i, alpha in enumerate(alphas):
            # Get data for this alpha
            theta = np.array(theta_dict[alpha])
            omega = np.array(omega_dict[alpha])
            
            # Wrap theta to [-pi, pi]
            theta = (theta + np.pi) % (2*np.pi) - np.pi
            
            # Convert to Manim points
            points = VMobject()
            manim_points = [axes.c2p(th, om) for th, om in zip(theta, omega)]
            points.set_points_as_corners(manim_points)
            
            # Style based on regime
            color = colors[i % len(colors)]
            if len(theta) < 10:  # Periodic
                points.set_stroke(color, 6, opacity=1)
                points.set_fill(color, opacity=0.8)
            else:  # Chaotic
                points.set_stroke(color, 2, opacity=0.9)
                points.set_fill(color, opacity=0.4)
            
            # Parameter label
            label = Text(f"$f = {alpha:.2f}$", font_size=28, color=color).to_corner(DR)
            
            # Animate appearance
            self.play(
                FadeIn(label),
                Create(points, run_time=2),
                run_time=2.5
            )
            
            all_points.add(points)
            self.wait(1)
            
            # Fade label but keep points
            self.play(FadeOut(label))
        
        # Final overview
        summary = Text("Chaos via Period-Doubling Cascade", font_size=32).to_edge(UP)
        self.play(
            all_points.animate.set_opacity(0.7),
            FadeIn(summary),
            run_time=2
        )
        self.wait(2)
        
        # Fade out for next phase
        self.play(
            FadeOut(all_points),
            FadeOut(summary),
            run_time=1
        )

class PoincareBuildup(Scene):
    """Bonus: Single chaotic attractor buildup animation."""
    
    def construct(self):
        # Show title
        title = Text("Chaotic Poincaré Section Buildup\n(f = 1.35)", font_size=32).to_edge(UP)
        self.play(Write(title))
        
        # Axes
        axes = Axes(
            x_range=[-np.pi, np.pi, np.pi/2],
            y_range=[-3, 3, 1],
            x_length=10,
            y_length=5
        ).shift(DOWN)
        
        self.play(Create(axes))
        
        # Load OR generate chaotic data
        # For demo, generate synthetic fractal-like data
        np.random.seed(42)
        n_points = 5000
        theta = np.random.normal(0, 1.2, n_points)
        omega = np.random.normal(0, 1.5, n_points)
        # Add fractal structure (simple)
        theta += 0.3 * np.sin(5*theta) + 0.1 * np.sin(12*theta)
        omega += 0.2 * np.cos(4*theta)
        
        theta = np.clip((theta + np.pi) % (2*np.pi) - np.pi, -np.pi, np.pi)
        
        # Progressive buildup
        points_group = VGroup()
        chunk_size = 100
        
        for i in range(0, n_points, chunk_size):
            chunk_theta = theta[i:i+chunk_size]
            chunk_omega = omega[i:i+chunk_size]
            
            chunk_points = VMobject()
            manim_pts = [axes.c2p(th, om) for th, om in zip(chunk_theta, chunk_omega)]
            chunk_points.set_points_as_corners(manim_pts)
            chunk_points.set_stroke(RED_C, 2, opacity=0.9)
            chunk_points.set_fill(RED_C, opacity=0.3)
            
            self.play(Create(chunk_points), run_time=0.5)
            points_group.add(chunk_points)
            
            if i % 500 == 0:
                self.add(Text(f"Points: {min(i+chunk_size, n_points)}", font_size=24).to_corner(DR))
        
        self.wait(2)
        self.play(FadeOut(points_group), FadeOut(title))

# Usage:
# 1. Run: python poincare_section_animation.py  (computes and stores data)
# 2. Run: manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\continuous_phase_portrait.py PhasePortraitScene -o my_animation.mp4