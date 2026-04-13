import numpy as np
from dataclasses import dataclass

# ============================================================
# 1. PARAMETERS AND EQUATION (Taylor form)
# ============================================================

@dataclass
class Params:
    beta: float = 0.25          # damping
    omega: float = 2/3          # drive frequency
    dt: float = 0.01            # time step


def rhs(t, y, p, gamma):
    """Taylor form: theta'' + beta*theta' + sin(theta) = gamma*cos(omega t)."""
    theta, omega = y
    dtheta = omega
    domega = -p.beta * omega - np.sin(theta) + gamma * np.cos(p.omega * t)
    return np.array([dtheta, domega])


# ============================================================
# 2. RK4 INTEGRATOR
# ============================================================

def rk4_step(f, t, y, dt, p, gamma):
    k1 = f(t,         y,           p, gamma)
    k2 = f(t+dt/2,    y+dt*k1/2,   p, gamma)
    k3 = f(t+dt/2,    y+dt*k2/2,   p, gamma)
    k4 = f(t+dt,      y+dt*k3,     p, gamma)
    return y + dt*(k1 + 2*k2 + 2*k3 + k4)/6


def integrate_steps(p, gamma, n_steps, y0, t0=0.0):
    """
    Integrate for exactly n_steps steps.
    Returns q(t), dq(t), and final state y_end.
    """
    q = np.empty(n_steps)
    dq = np.empty(n_steps)
    y = y0.copy()
    t = t0
    dt = p.dt

    for i in range(n_steps):
        q[i], dq[i] = y
        y = rk4_step(rhs, t, y, dt, p, gamma)
        t += dt

    return q, dq, y, t

# ============================================================
# 3. POINCARÉ SECTION
# ============================================================

def poincare_section(p, gamma, y0, t_transient=300.0, n_samples=200):
    """
    Compute Poincaré section by sampling once per drive period
    after a transient. This is phase-locked and matches the
    standard driven pendulum bifurcation diagrams.
    """
    T = 2.0 * np.pi / p.omega
    dt = p.dt
    steps_transient = int(t_transient / dt)
    steps_per_period = int(round(T / dt))

    # 1) Transient
    q_tr, dq_tr, y, t = integrate_steps(p, gamma, steps_transient, y0, t0=0.0)

    # 2) Poincaré sampling
    q_p = []
    dq_p = []

    for _ in range(n_samples):
        # integrate exactly one drive period
        q_tmp, dq_tmp, y, t = integrate_steps(p, gamma, steps_per_period, y, t0=t)
        theta, omega_v = y

        # wrap theta to [-pi, pi]
        theta = (theta + np.pi) % (2.0*np.pi) - np.pi

        q_p.append(theta)
        dq_p.append(omega_v)

    return np.array(q_p), np.array(dq_p), y, t



# ============================================================
# 4. LYAPUNOV EXPONENT
# ============================================================

def lyapunov(p, gamma, y0, t_transient=300, n_periods=800, delta0=1e-8):
    T = 2*np.pi / p.omega

    # transient
    q_tr, dq_tr = integrate_steps(p, gamma, t_transient, y0)
    y1 = np.array([q_tr[-1], dq_tr[-1]])
    y2 = y1 + np.array([delta0, 0.0])

    sum_log = 0.0
    t = 0.0
    steps_per_period = int(T / p.dt)
    total_steps = n_periods * steps_per_period

    for i in range(total_steps):
        y1 = rk4_step(rhs, t, y1, p.dt, p, gamma)
        y2 = rk4_step(rhs, t, y2, p.dt, p, gamma)
        t += p.dt

        if (i+1) % steps_per_period == 0:
            d = np.linalg.norm(y2 - y1)
            if d == 0:
                d = 1e-16
            sum_log += np.log(d / delta0)
            # renormalize
            y2 = y1 + (delta0 / d) * (y2 - y1)

    return sum_log / (n_periods * T)


# ============================================================
# 5. SWEEP OVER ALPHAS AND STORE EVERYTHING
# ============================================================

def compute_and_store(alphas, filename="C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\oscillatory_motion\\Driven_oscillation\\store_data.npz"):
    p = Params()
    data = {
        "q_full": {},
        "dq_full": {},
        "q_poincare": {},
        "dq_poincare": {},
        "lambda": {}
    }

    y0_base = np.array([0.0, 0.0])  # initial state

    for i, alpha in enumerate(alphas):
       
        progress = (i + 1) / len(alphas)
        bar_length = 12
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)

        print(rf"[{bar}]  {progress*100:5.1f}%   $\alpha$ = {alpha:.2f}", end="\r", flush=True)
            

      # 1) Full trajectory (for visualization, not too long)
        t_full = 150.0
        steps_full = int(t_full / p.dt)
        q_full, dq_full, _, _ = integrate_steps(p, alpha, steps_full, y0_base, t0=0.0)

        # 2) Poincaré section (phase-locked)
        q_p, dq_p, _, _ = poincare_section(
            p, alpha, y0_base, t_transient=300.0, n_samples=200
        )

        # store
        data["q_full"][alpha] = q_full
        data["dq_full"][alpha] = dq_full
        data["q_poincare"][alpha] = q_p
        data["dq_poincare"][alpha] = dq_p

    np.savez(filename, **data)
    print(f"\nSaved to {filename}")


# ============================================================
# 6. RUN
# ============================================================

if __name__ == "__main__":

    alphas = np.linspace(1.060, 1.087, 50)
    compute_and_store(alphas, "store_data.npz")
