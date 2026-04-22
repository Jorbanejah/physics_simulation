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
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)

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
        self.play(group.animate.to_corner(UR).shift(DOWN * 0.5))

    def create_phase_portrait_axes(self, trajectory):
        """Create 3D axes for trajectories and 2D axes for Poincaré sections."""
        #Proper range axis

        q_dict, dq_dict, t_dict = trajectory

        q_vals = np.concatenate(list(q_dict.values()))
        dq_vals = np.concatenate(list(dq_dict.values()))
        t_vals = np.concatenate(list(t_dict.values()))

        pad = 0.1

        q_min, q_max = q_vals.min() - pad, q_vals.max() + pad
        dq_min, dq_max = dq_vals.min() - pad, dq_vals.max() + pad
        t_min, t_max = t_vals.min() - pad, t_vals.max() + pad

        # 3D axes (affected by camera)
        axes3d = ThreeDAxes(
            x_range=[q_min, q_max, (q_max - q_min) / 6],
            y_range=[t_min, t_max, (t_max - t_min) / 6],
            z_range=[dq_min, dq_max, (dq_max - dq_min) / 6],
            x_length=5.5,
            y_length=5.5,
            z_length=5.5,
        ).shift(LEFT * 4)

        labels3d = axes3d.get_axis_labels(
            x_label="q",
            y_label="t_{mod T}",
            z_label=r"\dot{q}"
        )
        trajectories = Text("Trajectories", font_size=24)
    
        self.add(axes3d, labels3d, trajectories)

        # 2D axes (must stay flat)
        axes2d = Axes(
            x_range=[q_min, q_max, (q_max - q_min) / 6],
            y_range=[dq_min, dq_max, (dq_max - dq_min) / 6],
            x_length=5.5,
            y_length=5.5,
            axis_config={"include_tip": False, "font_size": 24},
            tips=False
        ).shift(RIGHT * 3)

        poincare = Text("Poincaré sections", font_size=24).next_to(axes2d, UP, buff=0.2)
        labels2d = axes2d.get_axis_labels(x_label="q", y_label=r"\dot{q}")

        # Keep 2D axes fixed in screen space
        self.add_fixed_in_frame_mobjects(axes2d, labels2d, poincare)

        self.play(
            FadeIn(axes3d), FadeIn(trajectories),
            FadeIn(axes2d), FadeIn(labels2d), FadeIn(poincare),
            run_time=1.5
        )
        self.wait(1)

        return axes3d, axes2d

    def show_phase_portraits(self, axes2d, axes3d, y_poin, y_traj, alphas):

        labels = [
            rf"1-Periodic: α = {alphas[0]:.3f}",
            rf"2-Periodic: α = {alphas[1]:.3f}",
            rf"Quasi-periodic: α = {alphas[2]:.3f}",
            rf"Chaotic: α = {alphas[3]:.3f}"
        ]

        all_curves = []
        all_points = []

        # 3D trajectories
        for alpha in alphas:
            curve = self.make_continuous_phase_plot(axes3d, y_traj, alpha)
            all_curves.append(curve)

        # 2D Poincaré points
        q_poin, dq_poin = y_poin

        for alpha in alphas:
            pts = [
                axes2d.c2p(q_poin[alpha][j], dq_poin[alpha][j])
                for j in range(len(q_poin[alpha]))
            ]
            dots = VGroup(*[Dot(pt, radius=0.05, color=YELLOW) for pt in pts])

            # Keep dots fixed in 2D
            self.add_fixed_in_frame_mobjects(dots)

            all_points.append(dots)

        # Animate each regime
        for i, (curve, points, label) in enumerate(zip(all_curves, all_points, labels)):

            text = Text(label, font_size=24).to_corner(DR).shift(LEFT * 1.5 + DOWN * 0.3) # cambiar al centro de la animacion abajo.
            self.add_fixed_in_frame_mobjects(text)

            self.play(
                Write(text),
                Create(curve, run_time=5),
                LaggedStart(*[FadeIn(dot, scale=1.5) for dot in points], lag_ratio=0.05),
                run_time=4
            )

            self.wait(1)

            self.play(FadeOut(text), FadeOut(curve), FadeOut(points))

        # Final 3D rotation
        self.begin_ambient_camera_rotation(rate=0.2)
        self.wait(3)
        self.stop_ambient_camera_rotation()

    def make_continuous_phase_plot(self, axes3d, y_traj, alpha):

        theta, omega, period = y_traj
        q = np.array(theta[alpha])
        dq = np.array(omega[alpha])
        T = np.array(period[alpha])

        points = [axes3d.c2p(q[i], dq[i], T[i]) for i in range(len(q))]

        curve = VMobject()  
        curve.set_points_smoothly(points)
        curve.set_stroke(RED_C, 4, opacity=0.9)

        return curve

# manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\traject_poincare_anim.py PhasePortraitScene -o my_animation.mp4