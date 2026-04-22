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
        axes3d, axes2d = self.create_phase_portrait_axes()

        # 3. Animate everything
        self.show_phase_portraits(axes2d, axes3d, y_poin, y_traj, alphas)

    def show_equation(self):
        """Display the driven pendulum equation in 2D."""
        title = Text("Driven pendulum equation").to_corner(UL)
        eq = MathTex(r"\ddot{q} + \beta \dot{q} + \sin(q) = \alpha \cos(\omega t)").scale(0.9)

        # Keep equation fixed in 2D
        self.add_fixed_in_frame_mobjects(title, eq)

        self.play(FadeIn(title, shift=UR))
        self.play(Write(eq))
        self.wait()
        self.play(FadeOut(title), FadeOut(eq))

    def create_phase_portrait_axes(self):
        """Create 3D axes for trajectories and 2D axes for Poincaré sections."""

        # 3D axes (affected by camera)
        axes3d = ThreeDAxes(
            x_range=(-np.pi, np.pi, 0.5),
            y_range=(-3, 3, 1),
            z_range=(-2, 2, 1),
            x_length=5.5,
            y_length=5.5,
            z_length=5.5,
        ).shift(LEFT * 3)

        trajectories = Text("Trajectories", font_size=24).next_to(axes3d, UP, buff=0.2)
        # 3D label rotates with camera
        self.add(axes3d, trajectories)

        # 2D axes (must stay flat)
        axes2d = Axes(
            x_range=[-np.pi, np.pi, np.pi/2],
            y_range=[-3, 3, 0.5],
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
            rf"1-Periodic: $\alpha = {alphas[0]:.3f}$",
            rf"2-Periodic: $\alpha = {alphas[1]:.3f}$",
            rf"Quasi-periodic: $\alpha = {alphas[2]:.3f}$",
            rf"Chaotic: $\alpha = {alphas[3]:.3f}$"
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
                for j in range(min(100, len(q_poin[alpha])))
            ]
            dots = VGroup(*[Dot(pt, radius=0.05, color=YELLOW) for pt in pts])

            # Keep dots fixed in 2D
            self.add_fixed_in_frame_mobjects(dots)

            all_points.append(dots)

        # Animate each regime
        for i, (curve, points, label) in enumerate(zip(all_curves, all_points, labels)):

            text = Tex(label).to_corner(UR).shift(LEFT * 2 + DOWN * i * 0.3)
            self.add_fixed_in_frame_mobjects(text)

            # Rotate only the 3D camera
            self.move_camera(
                phi=45 * DEGREES + i * 5 * DEGREES,
                theta=45 * DEGREES + i * 10 * DEGREES,
                run_time=0.5
            )

            self.play(
                Write(text),
                Create(curve, run_time=4),
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
        curve.set_stroke(BLUE_C, 4, opacity=0.9)

        return curve

# manim -pqh Desktop\Programas\Python\physics_simulation\oscillatory_motion\Driven_oscillation\traject_poincare_anim.py PhasePortraitScene -o my_animation.mp4