from manim import *
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
        # 3. Create poincare sections
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

        #Slide the equation away
        self.play(
            eq.animate.to_corner(UR)
        )

        self.new_title = Text("Poincare sections").to_corner(UL)
        self.play(Transform(title, self.new_title))

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
            x_range= axes.x_range,
            y_range= axes.y_range,
            x_length=8,
            y_length=5.5,
            background_line_style={
                "stroke_opacity": 0.12,
                "stroke_width": 1,
                "stroke_color": WHITE
            }
        ).move_to(axes)

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
        dots.set_opacity(1)
        self.add(dots) #ensure dots are in the scene

        # Static label object
        alpha_label = MathTex("").scale(1.2).to_corner(DR)
        behavior_label = Text("", font_size=24).next_to(alpha_label, DOWN)

        self.add(alpha_label, behavior_label)

        # Set initial label values
        alpha_label.become(MathTex(f"\\alpha = {alpha0:.2f}").scale(1.2).to_corner(DR))
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
            new_dots.set_opacity(0)
            self.add(new_dots)

            # Behavior classification
            if alpha < 1.065:
                behavior = "Normal"
            elif alpha < 1.070:
                behavior = "Erratic"
            elif alpha < 1.075:
                behavior = "Chaotic"
            else:
                behavior = "Strongly Chaotic"

            # Update labels IN PLACE
          
            new_alpha_label = MathTex(f"\\alpha = {alpha:.2f}").scale(1.2).to_corner(DR)
            new_behavior_label = Text(behavior, font_size=24).next_to(new_alpha_label, DOWN)

            # Animate smooth transitions
            self.play(
                dots.animate.set_opacity(0),
                new_dots.animate.set_opacity(1),
                Transform(alpha_label, new_alpha_label),
                Transform(behavior_label, new_behavior_label),
                rate_func=smooth,
                run_time=1.5
            )
            self.remove(dots)
            dots = new_dots
            self.wait(0.3)

        self.wait(2)

    def make_dots(self, axes, rk4_data, alpha):
        """Create dot group for Poincaré section at given alpha."""
        q = np.array(rk4_data["q"][alpha])
        dq = np.array(rk4_data["dq"][alpha])
        
        # Limit points for performance and clarity
        q = (q + np.pi) % (2*np.pi) - np.pi

        #Sorted points
        order = np.argsort(q)
        q = q[order]
        dq = dq[order]

        n_points = min(800, len(q))
        indices = np.linspace(0, len(q)-1, n_points, dtype=int)

        if alpha < 3:
            dots = VGroup(*[
                Dot(
                    axes.c2p(q[i], dq[i]), 
                    radius=0.04, 
                    color=BLUE_C
                )
                for i in indices
            ])
        elif alpha < 6:
            dots = VGroup(*[
                Dot(
                    axes.c2p(q[i], dq[i]), 
                    radius=0.04, 
                    color=YELLOW_C
                )
                for i in indices
            ])
        elif alpha < 12:
            dots = VGroup(*[
                Dot(
                    axes.c2p(q[i], dq[i]), 
                    radius=0.04, 
                    color=GREEN_C
                )
                for i in indices
            ])

        else:
            dots = VGroup(*[
                Dot(
                    axes.c2p(q[i], dq[i]), 
                    radius=0.04, 
                    color=RED_C
                )
                for i in indices
            ])
        
        return dots

#Run: manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\poincare_sections_anim.py PoincareSectionsScene -o my_animation.mp4