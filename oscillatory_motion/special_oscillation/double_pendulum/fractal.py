"""
This computes the "flip time" for each initial condition (θ₁, θ₂)
with both pendulums starting from rest (ω₁ = ω₂ = 0).

Reference: FIG. 2 shows the outcome regions for the double pendulum
where angles range from -π to π.

Color scheme:

Green: flips within 10 units

Red: 10-100 units

Purple: 100-1000 units

Blue: 1000-10000 units

White: doesn't flip within 10000 units

Black curve: energetically impossible (3cos(θ₁) + cos(θ₂) = 2)
"""

import numpy as np
import numba as nb
from numba import njit, prange
import matplotlib.pyplot as plt

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================

g = 9.81
m1 = 1.0
m2 = 1.0
L1 = 1.0
L2 = 1.0

flip_threshold = np.pi
t_max = 10000.0
dt = 0.01

# ============================================================
# ENERGY-BASED FLIP POSSIBILITY
# ============================================================

@njit
def can_flip(theta1, theta2):
    return 3*np.cos(theta1) + np.cos(theta2) >= 2.0

# ============================================================
# DOUBLE PENDULUM EQUATIONS (FAST)
# ============================================================

@njit
def derivatives(state):
    theta1, theta2, w1, w2 = state

    delta = theta1 - theta2
    s = np.sin(delta)
    c = np.cos(delta)

    M11 = (m1 + m2) * L1 * L1
    M12 = m2 * L1 * L2 * c
    M22 = m2 * L2 * L2

    F1 = -(m1 + m2) * g * L1 * np.sin(theta1) - m2 * L1 * L2 * w2 * w2 * s
    F2 = m2 * L1 * L2 * w1 * w1 * s - m2 * g * L2 * np.sin(theta2)

    det = M11 * M22 - M12 * M12

    a1 = (M22 * F1 - M12 * F2) / det
    a2 = (-M12 * F1 + M11 * F2) / det

    return np.array([w1, w2, a1, a2])

# ============================================================
# RK4 INTEGRATOR (JIT-COMPILED)
# ============================================================

@njit
def rk4_step(state, dt):
    k1 = derivatives(state)
    k2 = derivatives(state + 0.5 * dt * k1)
    k3 = derivatives(state + 0.5 * dt * k2)
    k4 = derivatives(state + dt * k3)
    return state + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

# ============================================================
# FLIP TIME SIMULATION (FAST)
# ============================================================

@njit
def simulate_flip(theta1_0, theta2_0):
    if not can_flip(theta1_0, theta2_0):
        return -2.0  # energetically forbidden

    state = np.array([theta1_0, theta2_0, 0.0, 0.0])
    t = 0.0

    for _ in range(int(t_max / dt)):
        state = rk4_step(state, dt)
        t += dt

        # Flip detection
        if abs(state[0]) >= flip_threshold or abs(state[1]) >= flip_threshold:
            return t

    return -1.0  # no flip

# ============================================================
# PARALLEL FRACTAL COMPUTATION
# ============================================================

@njit(parallel=True)
def compute_fractal_fast(N):
    theta_vals = np.linspace(-np.pi, np.pi, N)
    result = np.zeros((N, N))

    for i in prange(N):
        for j in range(N):
            result[j, i] = simulate_flip(theta_vals[i], theta_vals[j])

    return theta_vals, theta_vals, result

# ============================================================
# PLOT
# ============================================================

def plot_fractal(theta1, theta2, flip_times):
    vis = flip_times.copy()
    vis[vis <= 0] = np.nan
    log_vis = np.log10(np.clip(vis, 1, 1e10))

    plt.figure(figsize=(10, 8))
    plt.imshow(log_vis, origin='lower',
               extent=[-np.pi, np.pi, -np.pi, np.pi],
               cmap='viridis')
    plt.colorbar(label='log10(flip time)')
    plt.xlabel("θ₁")
    plt.ylabel("θ₂")
    plt.title("Fast Double Pendulum Flip-Time Fractal")
    plt.show()

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    N = 200  # 

    print("=" * 60)
    print("DOUBLE PENDULUM FRACTAL COMPUTATION")
    print("=" * 60)
    print(f"Grid size: {N} x {N}")
    print(f"Time limit: {t_max}")

    theta1, theta2, flips = compute_fractal_fast(N)
    plot_fractal(theta1, theta2, flips)