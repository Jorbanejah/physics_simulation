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

        # Title
        title = Text("Driven pendulum equation").to_corner(UL)

        # Main equation (dimensional form)
        eq = MathTex(
            r"\ddot{q} + \frac{\gamma}{m L^2 \frac{g}{L}} \dot{q} +  \sin(q) = \frac{F_0}{m L^2 \frac{g}{L}} \cos(\omega t)"
        ).scale(0.9)

        self.play(FadeIn(title, shift=DOWN))
        self.play(Write(eq))
        self.wait()

        # Nondimensional form   
        transform_eq = MathTex(
            r"\ddot{q} + \beta \dot{q} + \sin(q) = \alpha \cos(\omega t)"
        ).scale(0.9)

        # Make sure the new equation starts exactly where the old one is
        transform_eq.move_to(eq.get_center())

        self.play(Transform(eq, transform_eq))
        self.wait()

        # Slide the equation to the upper-right corner
        self.play(eq.animate.to_corner(UR))

        # Replace the title
        new_title = Text("Poincare sections").to_corner(UL)
        self.play(Transform(title, new_title))

        self.wait()

        
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
        

        self.play(
            Create(grid),
            Create(axes),
            FadeIn(x_label),
            FadeIn(y_label),
            Create(section_line),
            run_time=2
        )
        self.wait(0.5)

        return axes
    
    def show_poincare_sections(self, axes, y, alphas):
        """Animate 50 Poincaré sections in ~1.5 minutes, fading each out cleanly."""
    
        theta_dict, omega_dict = y
        colors = [BLUE_C, GREEN_C, ORANGE, RED_C, PURPLE, TEAL]

        total_time = 60  # 1 minutes
        time_per_alpha = total_time / len(alphas)

        # Create persistent alpha label
        alpha_label = Text(r"$\alpha$ = ", font_size=24)
        alpha_value = DecimalNumber(0, num_decimal_places=4, font_size=24)
        alpha_group = VGroup(alpha_label, alpha_value).arrange(RIGHT).to_corner(DR)
        self.add(alpha_group)
        
        # Filtering threshold
        eps_q = 1e-3
        eps_dq = 1e-3

        for i, alpha in enumerate(alphas):

            theta = np.array(theta_dict[alpha])
            omega = np.array(omega_dict[alpha])

            # === FILTER: remove points too close to zero ===
            if (np.abs(theta) < eps_q).all() and (np.abs(omega) < eps_dq).all():
                continue
            
            dots = VGroup(*[
            Dot(
                axes.c2p(th, om),
                radius=0.03,
                color=colors[i % len(colors)],
                stroke_width=0,
                fill_opacity=0.9
            )
                for th, om in zip(theta, omega)
            ])

            # Label for this alpha
            self.play(alpha_value.animate.set_value(alpha), run_time=0.3)

            # Show label + dots
            self.play(
                LaggedStart(*[FadeIn(d, scale=0.5) for d in dots],
                            lag_ratio=0.01,
                            run_time=time_per_alpha * 0.7),
                run_time=time_per_alpha,
                rate_func=linear
            )

            # Fade both out before next alpha
            self.play(
                FadeOut(dots),
                run_time=0.2
            )

        self.wait(1)


class PoincareBuildup(Scene):
    """Bonus: Single chaotic attractor buildup animation."""
    
    def construct(self):
        # Show title
        title = Text(r"Chaotic Poincaré Section Buildup\n($\alpha$ = 1.35)", font_size=32).to_edge(UP)
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
        # Generate synthetic fractal-like data
        np.random.seed(42)
        n_points = 100
        theta = np.random.normal(0, 1.2, n_points)
        omega = np.random.normal(0, 1.5, n_points)
        # Add fractal structure (simple)
        theta += 0.3 * np.sin(5*theta) + 0.1 * np.sin(12*theta)
        omega += 0.2 * np.cos(4*theta)
        
        theta = np.clip((theta + np.pi) % (2*np.pi) - np.pi, -np.pi, np.pi)
        
        # Progressive buildup
        chunk_size = 100
        
        for i in range(0, n_points, chunk_size):
            chunk_theta = theta[i:i+chunk_size]
            chunk_omega = omega[i:i+chunk_size]
            
            chunk_points = VMobject()
            manim_pts = [axes.c2p(th, om) for th, om in zip(chunk_theta, chunk_omega)]
            chunk_points.set_points_as_corners(manim_pts)
            chunk_points.set_stroke(RED_C, 2, opacity=0.9)
            
            self.play(Create(chunk_points), 
                run_time=0.3
            )

            self.play(FadeOut(chunk_points))
         
            if i % 500 == 0:
                self.add(Text(f"Points: {min(i+chunk_size, n_points)}", font_size=24).to_corner(DR))
        
        self.wait(2)

# Usage:
# 1. Run: python poincare_section_animation.py  (computes and stores data)
# 2. Run: manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\continuous_phase_portrait.py PhasePortraitScene -o my_animation.mp4