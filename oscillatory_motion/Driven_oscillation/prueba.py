
import numpy as np
from manim import *
from collections import defaultdict


class LyapunovBifurcation(Scene):

    def construct(self):
    
        data = np.load("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\data_bifurcation.npz", allow_pickle=True)

        #Load bifurcation
        q_p = data["bifur_q"] # length = 7500
        alphas_bi = data["alphas"] # length =7500

        #Load Lyapunov
        Lyapunov = data["Lyapunov"] # length 150
        alphas_L = np.linspace(1.060, 1.087, 150)
       
        bifur_dict = defaultdict(list)

        for A, q in zip(alphas_bi, q_p):
            bifur_dict[A].append(q)

        unique_alphas = sorted(bifur_dict.keys())

        #1. Show equation
        self.show_equation()

        #2. Create axes, labels, and grid
        ax_bif, ax_lyap = self.axes(alphas_L, q_p)

        #3. Animate Lyapunov and bifurcation diagram
        self.Bifurcation_Lyapunov(alphas_L, Lyapunov, bifur_dict, unique_alphas, ax_bif, ax_lyap)

    def show_equation(self):
        """Display the driven pendulum equation."""
        title = Text("Driven pendulum equation").to_corner(UL)
        eq = MathTex(r"\ddot{q} + \beta \dot{q} + \sin(q) = \alpha \cos(\omega t)").scale(0.9)

        self.play(FadeIn(title, shift=UR))
        self.play(Write(eq))
        self.wait()
        self.play(FadeOut(title), FadeOut(eq))

    def axes(self, alphas, q_p):
        """Create properly proportioned axes that fit the animation"""
        # BIFURCATION AXES - TOP
        ax_bif = Axes(
            x_range=[alphas[0], alphas[-1], 0.001],
            y_range=[min(q_p), max(q_p), 0.15],
            x_length=8.0,
            y_length=2.8,
            tips=False
        )

        bif_label = Text("Bifurcation Diagram (α)", font_size=24).next_to(ax_bif, UP, buff=0.2)

        # LYAPUNOV AXES - BOTTOM
        ax_lyap = Axes(
            x_range=[alphas[0], alphas[-1], 0.001],
            y_range=[-0.25, 0.15, 0.025],
            x_length=8.0,
            y_length=2.8,
            tips=False
        )

        lyap_label = Text("Lyapunov Exponent (λ)", font_size=24).next_to(ax_lyap, UP, buff=0.2)

        # STACK BOTH AXES VERTICALLY
        axes_group = VGroup(ax_bif, ax_lyap).arrange(DOWN, buff=0.3)
        labels_group = VGroup(bif_label, lyap_label).arrange(DOWN, buff=0.3)
        
        main_group = VGroup(labels_group, axes_group).move_to(ORIGIN)
        
        self.play(FadeIn(main_group))
        return ax_bif, ax_lyap

    def Bifurcation_Lyapunov(self, alphas_L, Lyapunov, bifur_dict, unique_alphas, ax_bif, ax_lyap):
        """Main animation loop with points, proper labels, and bifurcation lines"""

        min_len = min(len(unique_alphas), len(alphas_L))
        unique_alphas = unique_alphas[:min_len]
        alphas_L_loop = alphas_L[:min_len]

        # Create Lyapunov curve
        lyap_curve = VMobject().set_stroke(RED, 4)
        lyap_points = []

        # alpha label next to bifurcation plot (top)
        alpha_label = Text("α =", font_size=28).set_color(BLUE)
        alpha_value = DecimalNumber(
            alphas_L[0], num_decimal_places=4, font_size=28
        ).set_color(BLUE)
        alpha_group = VGroup(alpha_label, alpha_value).arrange(RIGHT, buff=0.1)
        alpha_group.next_to(ax_bif, RIGHT, buff=0.3).shift(UP*0.3)

        # lambda label next to Lyapunov plot (bottom)
        lambda_label = Text("λ =", font_size=28).set_color(RED)
        lambda_value = DecimalNumber(
            0, num_decimal_places=4, font_size=28
        ).set_color(RED)
        lambda_group = VGroup(lambda_label, lambda_value).arrange(RIGHT, buff=0.1)
        lambda_group.next_to(ax_lyap, RIGHT, buff=0.3)

        self.play(FadeIn(alpha_group), FadeIn(lambda_group))

        # Bifurcation detection
        bifurcation_alphas = self.detect_bifurcations(alphas_L, Lyapunov)
        print(f"Detected bifurcations at: {[f'{a:.4f}' for a in bifurcation_alphas]}")

        # Create bifurcation lines upfront (invisible)
        vlines_bif = VGroup()
        vlines_lyap = VGroup()
        for alpha_bif in bifurcation_alphas:
            vline_bif = ax_bif.get_vertical_line(
                ax_bif.c2p(alpha_bif, 0), 
                color=YELLOW, 
                stroke_width=6
            ).set_opacity(0)
            vline_lyap = ax_lyap.get_vertical_line(
                ax_lyap.c2p(alpha_bif, 0), 
                color=YELLOW, 
                stroke_width=6
            ).set_opacity(0)
            vlines_bif.add(vline_bif)
            vlines_lyap.add(vline_lyap)

        # MAIN ANIMATION LOOP
        current_bif_points = VGroup()
        
        for i, (alpha_b, alpha_L) in enumerate(zip(unique_alphas, alphas_L_loop)):
            # Update labels
            self.play(
                alpha_value.animate.set_value(alpha_L),
                lambda_value.animate.set_value(Lyapunov[i]),
                run_time=0.1
            )

            # Add Lyapunov point and update curve
            lyap_point = ax_lyap.c2p(alpha_L, Lyapunov[i])
            lyap_points.append(lyap_point)
            
            if len(lyap_points) > 1:
                lyap_curve.set_points_smoothly(lyap_points)

            # Create new bifurcation points for current alpha
            theta_vals = np.array(bifur_dict[alpha_b])
            new_points = VGroup(*[
                Dot(ax_bif.c2p(alpha_b, th), radius=0.012, color=BLUE)
                for th in theta_vals
            ])
            
            # Animate new points appearing
            self.play(
                FadeIn(new_points, scale=1.5),
                lyap_curve.animate.set_stroke(RED, 4),
                run_time=0.3
            )
            
            current_bif_points.add(new_points)

            # Show bifurcation lines when we reach those alphas
            for j, vline_bif in enumerate(vlines_bif):
                if abs(alpha_L - bifurcation_alphas[j]) < 0.0005:
                    self.play(
                        FadeIn(vline_bif),
                        FadeIn(vlines_lyap[j]),
                        run_time=0.4
                    )

        self.wait(3)

    def detect_bifurcations(self, alphas, lyapunov):
        """Detect bifurcation points: 1 to 2, 2 to 4, 4 to chaos"""
        bifurcations = []
        prev_lam = None
        prev_prev_lam = None
        counter = 0

        for i, lam in enumerate(lyapunov):
            alpha = alphas[i]
            
            # Detect 1→2 (local minimum)
            if prev_prev_lam is not None and counter == 0:
                if prev_lam < prev_prev_lam and prev_lam < lam:
                    bifurcations.append(alpha)
                    counter += 1
            
            # Detect period doubling and chaos (sign changes)
            if prev_lam is not None and counter < 3:
                if prev_lam < 0 and lam > 0:
                    bifurcations.append(alpha)
                    counter += 1
            
            prev_prev_lam = prev_lam
            prev_lam = lam
        
        return bifurcations[:3]  # Return first 3 bifurcations