import numpy as np
import matplotlib.pyplot as plt

data = np.load("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\poincare_sections.npz", allow_pickle=True)

rk4_data = data["rk4"].item()

alphas = sorted(rk4_data["q_steady"].keys())

all_alphas = []
all_q = []
for alpha in alphas:
    q = np.array(rk4_data["q_steady"][alpha])

    # -pi to pi
    q = (q + np.pi) % (2*np.pi) - np.pi

    order = np.argsort(q)
    all_alphas.extend(len(q) * [alpha])
    all_q.extend(q)

plt.figure(figsize=(10, 6))
plt.scatter(x = all_alphas, y = all_q, s = 5)
plt.title("Bifurcation diagram")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\theta$")
plt.show()
