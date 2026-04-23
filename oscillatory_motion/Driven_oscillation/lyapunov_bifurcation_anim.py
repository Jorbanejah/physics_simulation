import numpy as np
from manim import *
from collections import defaultdict
#This program takes ~5/7 minutes to run. It animates 300

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

        # --- VALUE TRACKERS ---
        alpha_tracker = ValueTracker(alphas_L[0])
        lambda_tracker = ValueTracker(Lyapunov[0])

        # --- LABELS ---
           
        alpha_label = MathTex(r"\alpha =", font_size=24)
        alpha_label.next_to(ax_bif, RIGHT, buff=0.5)

        alpha_value = always_redraw(
            lambda: DecimalNumber(
                alpha_tracker.get_value(),
                num_decimal_places=4,
                font_size=24
            ).next_to(alpha_label, RIGHT, buff=0.2)
        )

        self.add(alpha_label, alpha_value)

        lambda_label = MathTex(r"\lambda =", font_size=24)  
        lambda_label.next_to(ax_lyap, RIGHT, buff=0.5)

        lambda_value = always_redraw(
            lambda: DecimalNumber(
                lambda_tracker.get_value(),
                num_decimal_places=4,
                font_size=24
            ).next_to(lambda_label, RIGHT, buff=0.2)
        )

        self.add(lambda_label, lambda_value)

        # --- LYAPUNOV CURVE (grows dynamically) ---
        lyap_curve = VMobject().set_stroke(RED_C, 3)

        def update_lyap_curve(mob):
            # Compute all points up to current alpha
            current_alpha = alpha_tracker.get_value()
            mask = alphas_L <= current_alpha
            
            pts = [ax_lyap.c2p(a, l) for a, l in zip(alphas_L[mask], Lyapunov[mask])]
            if len(pts) > 1:
                mob.set_points_as_corners(pts)

        lyap_curve.add_updater(update_lyap_curve)
        self.add(lyap_curve)

        # --- BIFURCATION POINTS (only show those with alpha <= current) ---
        all_bif_points = VGroup()

        for a in unique_alphas:
            for th in bifur_dict[a]:
                dot = Dot(ax_bif.c2p(a, th), radius=0.015, color=BLUE_C)
                dot.alpha_value = a  # store alpha for filtering
                all_bif_points.add(dot)

        def update_bif_points(group):
            current_alpha = alpha_tracker.get_value()
            for dot in group:
                dot.set_opacity(1 if dot.alpha_value <= current_alpha else 0)

        all_bif_points.add_updater(update_bif_points)
        self.add(all_bif_points)

        # --- LAMBDA TRACKER UPDATER ---
        dummy = Dot(radius=0.0001, color=WHITE)

        def update_lambda(mob):
            current_alpha = alpha_tracker.get_value()
            idx = np.searchsorted(alphas_L, current_alpha)
            idx = np.clip(idx, 0, len(Lyapunov)-1)
            lambda_tracker.set_value(Lyapunov[idx])

        dummy.add_updater(update_lambda)
        self.add(dummy)

        # --- MAIN ANIMATION ---
        self.play(
           alpha_tracker.animate.set_value(alphas_L[-1]),
            run_time=25,
            rate_func=linear
        )

        self.wait(2)
    
        
# manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\lyapunov_bifurcation_anim.py LyapunovBifurcation -o my_animation.mp4