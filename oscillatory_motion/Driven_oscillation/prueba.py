from manim import *
import numpy as np

class PhasePortraitAxesDemo(ThreeDScene):
    """Standalone demo showing corrected 3D + 2D axes layout."""

    def construct(self):

        # Fake data just to generate ranges
        q_vals = np.linspace(-2, 2, 100)
        dq_vals = np.linspace(-3, 3, 100)
        t_vals = np.linspace(0, 1, 100)

        trajectory = (
            {"demo": q_vals},
            {"demo": dq_vals},
            {"demo": t_vals},
        )
        self.set_camera_orientation(phi=65 * DEGREES, theta=60 * DEGREES)
        
        axes3d, axes2d = self.create_phase_portrait_axes(trajectory)

        # Rotate camera to show 3D axes move while 2D axes stay fixed
       

    def create_phase_portrait_axes(self, trajectory):
        """Corrected version: 3D axes in world space, 2D axes fixed in screen space."""

        q_dict, dq_dict, t_dict = trajectory

        q_vals = np.concatenate(list(q_dict.values()))
        dq_vals = np.concatenate(list(dq_dict.values()))
        t_vals = np.concatenate(list(t_dict.values()))

        # Scale t_mod_T for visibility
        t_scale = 5
        t_vals_scaled = t_vals * t_scale

        pad = 0.1

        q_min, q_max = q_vals.min() - pad, q_vals.max() + pad
        dq_min, dq_max = dq_vals.min() - pad, dq_vals.max() + pad
        t_min, t_max = t_vals_scaled.min() - pad, t_vals_scaled.max() + pad

        # --------------------------
        # 3D AXES (world space)
        # --------------------------
        axes3d = ThreeDAxes(
            x_range=[q_min, q_max, (q_max - q_min) / 6],
            y_range=[t_min, t_max, (t_max - t_min) / 6],
            z_range=[dq_min, dq_max, (dq_max - dq_min) / 6],
            x_length=5.5,
            y_length=5.5,
            z_length=5.5,
        ).to_corner(RIGHT)

        labels3d = axes3d.get_axis_labels(
            x_label="q",
            y_label="t_{mod T}",
            z_label=r"\dot{q}"
        )

        self.add(axes3d, labels3d)

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
        title2d = Text("Poincaré sections", font_size=24).next_to(axes2d, UP, buff=0.2)

        self.add_fixed_in_frame_mobjects(axes2d, labels2d, title2d)

        # --------------------------
        # Animate appearance
        # --------------------------
        self.play(
            FadeIn(axes3d), FadeIn(labels3d),
            FadeIn(axes2d), FadeIn(labels2d), FadeIn(title2d),
            run_time=1.5
        )

        return axes3d, axes2d
