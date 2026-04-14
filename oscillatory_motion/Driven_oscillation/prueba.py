import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
beta = 0.5
w = 2/3
T = 2 * np.pi / w

t_transient = 100
t_steady = 50

# Correct equation: solve_ivp expects f(t, y)
def driven_equation(t, y, A):
    theta, omega = y
    return [omega, -beta * omega - np.sin(theta) + A * np.cos(w * t)]

# Initial condition
y0 = np.array([0.0, 0.0])

As = np.linspace(1.060, 1.087, 150)

results_A = []
results_theta = []

for A in As:

    # Integrate full transient + steady window
    sol = solve_ivp(
        driven_equation,
        [0, (t_transient + t_steady) * T],
        y0,
        args=(A,),
        dense_output=True,   # IMPORTANT: allows exact sampling at multiples of T
        max_step=0.05        # keeps integration stable
    )

    # Sample exactly at multiples of T
    t_eval = np.arange(t_transient * T, (t_transient + t_steady) * T, T)
    y_samples = sol.sol(t_eval)
    theta_samples = y_samples[0]

    # Wrap angle to [-π, π]
    theta_wrapped = (theta_samples + np.pi) % (2*np.pi) - np.pi

    # Store results
    results_A.extend([A] * len(theta_wrapped))
    results_theta.extend(theta_wrapped)

    # Update initial condition for next A (follows attractor branch)
    y0 = sol.y[:, -1]

# Plot
plt.figure(figsize=(10,6))
plt.scatter(results_A, results_theta, s=5, color="black")
plt.xlabel("A")
plt.ylabel("θ (wrapped)")
plt.title("Driven Pendulum Bifurcation Diagram")
plt.show()
