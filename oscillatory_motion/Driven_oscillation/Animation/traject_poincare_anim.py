from manim import *
import numpy as np

class PhasePortraitScene(ThreeDScene):
    """Main animation scene showing Poincaré section evolution."""

    def construct(self):

        # Load data
        data = np.load(
            "C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\data_poin_traj.npz",
            allow_pickle=True
        )

        q_poincare = data["q_poincare"].item()
        dq_poincare = data["dq_poincare"].item()

        q_trajectory = data["q_trajectory"].item()
        dq_trajectory = data["dq_trajectory"].item()
        period = data["t_mod_T"].item()

        y_poin = q_poincare, dq_poincare
        y_traj = q_trajectory, dq_trajectory, period

        alphas = sorted(dq_poincare.keys())

        # Set 3D camera orientation (only affects 3D objects)
        self.set_camera_orientation(phi=75 * DEGREES, theta=45 * DEGREES)

        # 1. Show equation (2D)
        self.show_equation()

        # 2. Create axes (3D + 2D)
        axes3d, axes2d = self.create_phase_portrait_axes(y_traj)

        # 3. Animate everything
        self.show_phase_portraits(axes2d, axes3d, y_poin, y_traj, alphas)

    def show_equation(self):
        """Display the driven pendulum equation in 2D."""
        title = Text("Driven pendulum equation").to_corner(UL)
        eq = MathTex(r"\ddot{q} + \beta \dot{q} + \sin(q) = \alpha \cos(\omega t)").scale(0.9)

        group = VGroup(title, eq).arrange(DOWN, aligned_edge=LEFT)
        group.to_corner(UL)
        # Keep equation fixed in 2D
        self.add_fixed_in_frame_mobjects(group)

        self.play(FadeIn(title, shift=DOWN))
        self.play(Write(eq))
        self.wait()
    
        #Slide the equation away
        self.play(FadeOut(eq, title))

        PS =Text("Poincare sections")
        Traj = Text("Trajectories")
        
        self.add_fixed_in_frame_mobjects(PS, Traj)
        PS.to_corner(DR)
        Traj.to_corner(DL)
        self.play(FadeIn(PS, Traj))
        
    def create_phase_portrait_axes(self, trajectory):
        """Corrected version: 3D axes in world space, 2D axes fixed in screen space."""

        q_dict, dq_dict, t_dict = trajectory

        q_vals = np.concatenate(list(q_dict.values()))
        dq_vals = np.concatenate(list(dq_dict.values()))
        t_vals = np.concatenate(list(t_dict.values()))

        pad = 0.1

        q_min, q_max = q_vals.min() - pad, q_vals.max() + pad
        dq_min, dq_max = dq_vals.min() - pad, dq_vals.max() + pad
        t_min, t_max = t_vals.min() - pad, t_vals.max() + pad

        # --------------------------
        # 3D AXES (world space)
        # --------------------------
        axes3d = ThreeDAxes(
            x_range=[-5, 5, (dq_max - dq_min) / 6],
            y_range=[-5, 5, (q_max - q_min)/ 6],
            z_range=[t_min, t_max, (t_max - t_min) / 6],
            x_length=7.5,
            y_length=7.5,
            z_length=3.5,
        ).to_corner(UR).shift(LEFT * 0.8)


        self.add(axes3d)

        # --------------------------
        # 2D AXES (fixed in frame)
        # --------------------------
        axes2d = Axes(
            x_range=[q_min, q_max, (q_max - q_min) / 6],
            y_range=[dq_min, dq_max, (dq_max - dq_min) / 6],
            x_length=5.5,
            y_length=5.5,
            axis_config={"include_tip": False, "font_size": 24},
            tips=False
        )

        axes2d.to_corner(UR).shift(LEFT * 0.5)

        labels2d = axes2d.get_axis_labels(x_label="q", y_label=r"\dot{q}")
        

        self.add_fixed_in_frame_mobjects(axes2d, labels2d)

        # --------------------------
        # Animate appearance
        # --------------------------
        self.play(
            FadeIn(axes3d),
            FadeIn(axes2d), FadeIn(labels2d),
            run_time=1.5
        )

        return axes3d, axes2d

    def show_phase_portraits(self, axes2d, axes3d, y_poin, y_traj, alphas):
        
        q_poin, dq_poin = y_poin
        q_traj, dq_traj, t_traj = y_traj

       # Dictionary to store everything together
        data = {}

    # ---------------------------------------------------------
    # 1. Precompute curves + poincaré points for each α
    # ---------------------------------------------------------
        for alpha in alphas:

            # --- 3D curve ---
            q = np.array(q_traj[alpha])
            dq = np.array(dq_traj[alpha])
            T = np.array(t_traj[alpha])

            points3d = [axes3d.c2p(q[i], dq[i], T[i]) for i in range(len(q))]
            curve = VMobject().set_points_smoothly(points3d)
            curve.set_stroke(RED_C, 4, opacity=0.9)

            # --- 2D poincaré points ---
            pts2d = []
            for qv, dqv in zip(q_poin[alpha], dq_poin[alpha]):
                pt = axes2d.c2p(qv, dqv)
                dot = Dot(pt, radius=0.05, color=YELLOW)
                pts2d.append(dot)

            data[alpha] = {
                "curve": curve,
                "points": pts2d
            }

    # ---------------------------------------------------------
    # 2. Animate each α regime
    # ---------------------------------------------------------
        labels = [
            rf"1-Periodic: α = {alphas[0]:.3f}",
            rf"2-Periodic: α = {alphas[1]:.3f}",
            rf"Quasi-periodic: α = {alphas[2]:.3f}",
            rf"Chaotic: α = {alphas[3]:.3f}"
        ]

        for alpha, label in zip(alphas, labels):

            curve = data[alpha]["curve"]
            poincare_points = data[alpha]["points"]

            # Label for this regime
            text = Text(label, font_size=24).to_corner(UL)
            self.add_fixed_in_frame_mobjects(text)

            # Animate curve + poincaré points
            self.play(
                Write(text),
                Create(curve, run_time=5),
            )

            # Point appearance
            points_group = VGroup(*poincare_points)
            self.add_fixed_in_frame_mobjects(points_group)
            self.play(FadeIn(points_group, scale=1.5), run_time=0.5)

            self.wait(1)

            # Fade out everything for next regime
            self.play(FadeOut(text), FadeOut(curve), FadeOut(points_group))


# manim -pqh traject_poincare_anim.py PhasePortraitScene -o my_animation.mp4