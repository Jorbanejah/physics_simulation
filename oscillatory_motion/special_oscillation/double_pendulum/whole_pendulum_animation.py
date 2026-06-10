import numpy as np
import matplotlib.pyplot as plt
import os
from manim import *
from double_pendulum import DoublePendulumSimulator
from typing import Dict, Tuple

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

class SimulateDoublePendulum(Scene):
    def contruct(self):

        trajectory1, trajectory2, trajectory3 = load_data()

        self.show_equation()

        #2. Create axes, labels and grid

        axes = self.create_axes()

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

        self.second_title = Tex("Animation of three double pendulum")
        self.play(Transform(title, self.second_title))

    def create_axes(self):
        pass

    def animation(self, axes, trajectories1, trajecories2, trajectories3):
        pass