import numpy as np
import matplotlib.pyplot as plt

data = np.load("C:\\Users\\JORGE\\store_data.npz", allow_pickle=True)
q_p = data["q_poincare"].item()

all_A = []
all_theta = []

for A in sorted(q_p.keys()):
    qvals = np.array(q_p[A])
    all_A.extend([A] * len(qvals))
    all_theta.extend(qvals)

plt.figure(figsize=(10,6))
plt.scatter(all_A, all_theta, s=5, color="black")
plt.xlabel("A")
plt.ylabel("θ (wrapped)")
plt.title("Driven Pendulum Bifurcation Diagram")
plt.show()
