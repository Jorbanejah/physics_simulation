import numpy as np
import matplotlib.pyplot as plt
import os
from manim import *
from typing import Dict, Tuple, Sequence

def load_data(filename:str = "trajectories.py")-> Tuple[Dict, Dict, Dict]:
    """
    Load precomputed simulation data and check whether the data
    """
    directory = os.getcwd()
    route = os.path.join(directory, filename)
    if not os.path.isfile(route):
        raise FileNotFoundError(f"Data file not found in {route}. Run compute_data.py first.")

    data = np.load(route, allow_pickle= True)

    trajectory1 = data["trajectory1"].item()
    trajectory2 = data["trajectory2"].item()
    trajectory3 = data["trajectory3"].item()

    return trajectory1, trajectory2, trajectory3

def x_y_range(trajectory1:Sequence[float], trajectory2: Sequence[float], trajectory3: Sequence[float])-> Tuple[float, float]:
            
            min_x = 0
            max_x = 0 
            min_y = 0
            max_y = 0
            for traj in [trajectory1, trajectory2, trajectory3]:

                menor_x = min(traj["x2"])
                menor_y = min(traj["y2"])
                mayor_x = max(traj["x2"])
                mayor_y = max(traj["y2"])

                if menor_x < min_x:
                    min_x = menor_x

                if menor_y < min_y:
                    min_y = menor_y

                if mayor_x < max_x:
                    max_x = mayor_x

                if mayor_y < max_y:
                    max_y = mayor_y

            return min_x, max_x, min_y, max_y

class SimulateDoublePendulum(Scene):
    def contruct(self):

        trajectory1, trajectory2, trajectory3 = load_data()

        self.show_equation()

        #2. Create axes, labels and grid
        
        range_axes = x_y_range(trajectory1= trajectory1 ,trajectory2=trajectory2, trajectory3=trajectory3)

        axes = self.create_axes(range_axes)

        #3. Create animation

        self.animation(axes = axes, trajectories1= trajectory1, trajecories2= trajectory2,  trajectories3=trajectory3)

    def show_equation(self):
        """
        Display the double pendulum equation
        """

        title = Text("Double pendulum equation").to_corner(UL)

        eq = MathTex(
            r"""
            \begin{aligned}
            M_{11} &= (m_1 + m_2) L_1^2, &
            M_{12} &= m_2 L_1 L_2 \cos(\theta_1 - \theta_2), &
            M_{21} &= m_2 L_1 L_2 \cos(\theta_1 - \theta_2), &
            M_{22} &= m_2 L_2^2 \\
            
            F_1 &= -(m_1 + m_2) g L_1 \sin\theta_1
            - m_2 L_1 L_2 \,\omega_2^2 \sin(\theta_1 - \theta_2) \\
            F_2 &= m_2 L_1 L_2 \,\omega_1^2 \sin(\theta_1 - \theta_2)
            - m_2 g L_2 \sin\theta_2
            \end{aligned}
            """
            )
        
        self.play(FadeIn(title, shift = DOWN))
        self.play(Write(eq2))
        self.wait()

        eq2 = MathTex(
            r"""
            \begin{aligned}
            \theta_1' &= \omega_1 \\
            \theta_2' &= \omega_2 \\
            \omega_1' &= \frac{M_{22} F_1 - M_{12} F_2}{M_{11} M_{22} - M_{12}^2} \\
            \omega_2' &= \frac{-M_{12} F_1 + M_{11} F_2}{M_{11} M_{22} - M_{12}^2}
            \end{aligned}
            """
        ).scale(0.9)

        self.wait()

        self.play(Transform(eq, eq2))

        self.second_title = Tex("Animation of three double pendulums")
        self.play(Transform(title, self.second_title))

    def create_axes(self, range_axes: Sequence[float]):
        "Display the axes for double pendulum in 2d"

        pad = 0.1
        axes2d = Axes(
            x_range = [range_axes[0] + pad, range_axes[1]+ pad, (range_axes[0] + range_axes[1]) / 6],
            y_range = [range_axes[2]+ pad, range_axes[3] + pad, (range_axes[2] + range_axes[3]) / 6],
            x_length=5.5,
            y_length= 5.5,
            axis_config={"index_tip": False, "font_size": 24},
            tips = False
        )

        labels2d = axes2d.get_axis_labels(x_label="x", y_label= "y")

        self.play(
             FadeIn(axes2d),
             FadeIn(labels2d),
             run_time = 1.5
        )
    
        return axes2d
    
    def animation(self, axes, trajectories1, trajecories2, trajectories3):
        
        #Precompute trajectories

        for traj in (trajectories1, trajecories2, trajectories3):
             
            x1, y1, x2, y2 = traj["x1"], traj["y1"], traj["x2"], traj["y2"]

            points = [axes.c2p(x1[i], y1[i], x2[i], y2[i]) for i in range(len(x2))]

            curve = VMobject().set_points_smoothly(points=points)
            

