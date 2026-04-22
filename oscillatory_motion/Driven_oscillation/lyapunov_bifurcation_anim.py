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

        unique_alphas = sorted(bifur_dict.keys()) # Now we have 150 values

        #1. Show equation
        self.show_equation()

        #2. Create axes, labels, and grid
        ax_bif, ax_lyap = self.axes(alphas_L, q_p)

        #3. Animate Lyapunov and bifurcation diagram
        self.Bifurcation_Lyapunov(alphas_L, Lyapunov, bifur_dict, unique_alphas, ax_bif, ax_lyap)

    def show_equation(self):

        """Display the driven pendulum equation."""

        #Create title
        title = Text("Driven pendulum equation").to_corner(UL)

        #Main equation   
        eq =  MathTex(
            r"\ddot{q} + \beta \dot{q} + \sin(q) = \alpha \cos(\omega t)"
        ).scale(0.9)

        self.play(FadeIn(title, shift=UR))

        self.play(Write(eq))

        self.wait()

        self.play(FadeOut(title),FadeOut(eq))

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

        lyap_curve = VMobject().set_stroke(RED, 3)
        lyap_points = []
        self.add(lyap_curve)

        total_time = 25
        time_per_alpha = total_time / len(alphas_L)

        # Persistent alpha and lambda label
        alpha_label = Text(r"$\alpha =$", font_size=24)
        alpha_value = DecimalNumber(0, num_decimal_places=4, font_size=24)
        alpha_group = VGroup(alpha_label, alpha_value).arrange(RIGHT)
        alpha_group.next_to(ax_bif, RIGHT, buff=0.5)

        self.add(alpha_group)

        lambda_label = Text(r"$\lambda =$", font_size=24)
        lambda_value = DecimalNumber(0, num_decimal_places=4, font_size=24)
        lambda_group = VGroup(lambda_label, lambda_value).arrange(RIGHT)
        lambda_group.next_to(ax_lyap, RIGHT, buff=0.5)

        self.add(lambda_group)

        # Store previous Lyapunov exponent
        prev_lam = None
        prev_prev_lam = None
        counter = 0

        # --- MAIN LOOP ---
        for alpha_b, alpha_L in zip(unique_alphas, alphas_L):

            # Update α label
            self.play(alpha_value.animate.set_value(alpha_L), run_time=0.2)

            # --- BIFURCATION POINTS ---
            theta_vals = np.array(bifur_dict[alpha_b])
            bif_points = VGroup(*[
                Dot(ax_bif.c2p(alpha_b, th), radius=0.015, color=BLUE)
                for th in theta_vals
            ])

            # --- LYAPUNOV POINT ---
            lam = Lyapunov[np.where(alphas_L == alpha_L)[0][0]]
            lyap_points.append(ax_lyap.c2p(alpha_L, lam))
            lyap_curve.set_points_as_corners(lyap_points)

            self.play(lambda_value.animate.set_value(lam), run_time =0.2)

            triggered = False

            # 1) Detect 2→4 and 4→chaos (sign change)
            if prev_lam is not None:
                if prev_lam < 0 and lam > 0 and counter < 4:
                    triggered = True
                    counter +=1

            # 2) Detect 1→2 (local minimum)
            if prev_prev_lam is not None:
                if prev_lam < prev_prev_lam and prev_lam < lam and counter == 0:
                    # local minimum == first bifurcation
                    triggered = True
                    counter +=1

            if triggered:
                vline_bif = ax_bif.get_vertical_line(
                    ax_bif.c2p(alpha_L, 0), color=YELLOW
                ).set_stroke(width=2)

                vline_lyap = ax_lyap.get_vertical_line(
                    ax_lyap.c2p(alpha_L, 0), color=YELLOW
                ).set_stroke(width=2)

                self.play(
                    FadeIn(vline_bif),
                    FadeIn(vline_lyap),
                    run_time=0.4
                )

            prev_prev_lam = prev_lam
            prev_lam = lam

            # Animate both
            self.play(
                FadeIn(bif_points, scale=0.5),
                lyap_curve.animate.set_stroke(RED, 3),
                run_time=time_per_alpha,
                rate_func=linear
            )

        self.wait(2)
        
# manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\lyapunov_bifurcation_anim.py LyapunovBifurcation -o my_animation.mp4