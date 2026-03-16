import numpy as np
from Damping_vibration import Spring
import matplotlib.pyplot as plt

def draw(t, t_max, spr):
    """
    Simulates the motion of a damped spring system over time and calculates energies.

    Parameters:
        t (float): Initial time.
        t_max (float): Maximum simulation time.
        spr (Spring): An instance of the Spring class representing the system.

    Returns:
        tuple: Lists of potential energy (Ep), kinetic energy (Ek), work done by damping (Wp),
               total energy (Et), time points (times), positions (y_axes), and velocities (vy_axes).
    """

    y = spr.y
    v =spr.v
    dt = spr.dt

    Ep = []
    Ek = []
    Wp = []
    times = [t]
    y_axes = [y]
    vy_axes= [v]

    Ek0, Ep0, Wp0, _ = spr.energy(y, v, t)
    Ep.append(Ep0)
    Ek.append(Ek0)
    Wp.append(Wp0)

    relaxation_time = None
    found = False

    time_points = np.arange(t, t_max, dt)

    for t in time_points:

        y, vy= spr.rk4(dt)

        # Store the results

        y_axes.append(y)
        times.append(t)
        vy_axes.append(vy)

        if (not found) and abs(y) <= abs(y_axes[0]) * np.exp(-1):
            relaxation_time = t
            found = True

        #Energies
        Ek0, Ep0, Wp0, _ = spr.energy(y, vy, spr.dt)
        Ep.append(Ep0)
        Ek.append(Ek0)
        Wp.append(Wp0 + Wp[-1])

        spr.y = y
        spr.v = vy

    Et = np.array(Ep) + np.array(Ek) + np.array(Wp)

    return Ep, Ek, Wp, Et, times, y_axes, vy_axes, relaxation_time

k_under = np.linspace(0.5, 2.9, 10)
k = [3, 8, 18]
m = 2
gamma = [2, 8, 13]
y0 = -2
vy0 = 1
# [1, 2, 3] -----> omega
# [0.5, 2, >3] --------> beta
# beta < omega
# beta = omega
# beta> omega

results = {
    'k': k,
    'y': {},
    'vy': {},
    'Ep': {},
    'Ek': {},
    'Wp': {},
    'Et': {},
    'times': {}
}
for i, values in enumerate(k):
    t_max = 15
    t = 0
    spr = Spring(y0, vy0, m, values, gamma[i], dt = 0.01, t_max =15)

    Ep, Ek, Wp, Et, times, y_list, vy_list, _ = draw(t, t_max, spr)

    results['y'][values] = y_list
    results['vy'][values] = vy_list
    results['Ep'][values] = Ep
    results['Ek'][values] = Ek
    results['Wp'][values] = Wp
    results['Et'][values] = Et
    results['times'][values] = times

##
# --------- Graphics ---------
##


# ---------------------------------------------------------
# 1. Colormap
# ---------------------------------------------------------
cmap = plt.colormaps["viridis"]
colors = [cmap(i / len(k)) for i in range(len(k))]
plt.figure(figsize=(10, 6))

for i, values in enumerate(k):
    plt.plot(results['vy'][values],
             results['y'][values],
             color=colors[i],
             lw=2,
             label=rf"$k={values}$ N/m")
plt.xlim(-1, 1)
plt.ylim(-4, 0)
plt.axhline(spr.y0, color="black", lw=1, linestyle="--", label="Equilibrium")
plt.xlabel("Time $t$ (s)")
plt.ylabel(" $y$ (m)")
plt.title("Damping oscillator")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
#plt.savefig("\\colormap_spring.png", dpi=300, bbox_inches='tight')
plt.show()


# ---------------------------------------------------------
# 2. TRAJECTORY GRAPHICS
# ---------------------------------------------------------
plt.figure(figsize=(10, 6))

for i, values in enumerate(k):
    plt.plot(results['times'][values],
             results['y'][values],
             color=colors[i],
             label=rf"$\omega_0^2 = {values/m:.1f}\,\mathrm{{rad^2/s^2}}$")

plt.axhline(spr.y0, color="black", lw=0.8, linestyle="--", label="Equilibrium")
plt.xlabel("Time $t$ (s)")
plt.ylabel(" $y$ (m)")
plt.title("Damping oscillator")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
#plt.savefig("\\trajectory_spring.png", dpi=300, bbox_inches='tight')
plt.show()

# ---------------------------------------------------------
# 3. Energy graphics
# ---------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(18, 6), tight_layout=True)

for i, values in enumerate(k):
    ax = axes[i]
    ax.plot(results['times'][values], results['Ek'][values], label="Kinetic energy $E_k$")
    ax.plot(results['times'][values], results['Ep'][values], label="Potencial energy $E_p$")
    ax.plot(results['times'][values], results['Wp'][values], label= "Dissipated energy $W_p$")
    ax.plot(results['times'][values], results['Et'][values], label="Total energy $E_t$", lw=2)

    ax.set_xlabel("Time $t$ (s)")
    ax.set_ylabel("Energy (J)")
    ax.set_title(rf"Energy for $k = {values}$ N/m")
    ax.grid(alpha=0.3)
    ax.legend()

#plt.savefig("\\energy_spring.png", dpi=300, bbox_inches='tight')
plt.show()
# ---------------------------------------------------------
# 4. Numerical error
# ---------------------------------------------------------
def analytic_solution(t, y0, v0, beta, omega0):
    if beta < omega0:  
        # Under-damped
        omega_d = np.sqrt(omega0**2 - beta**2)
        A = y0
        B = (v0 + beta*y0) / omega_d
        return np.exp(-beta*t) * (A*np.cos(omega_d*t) + B*np.sin(omega_d*t))

    elif beta == omega0:
        # critically damped
        return np.exp(-beta*t) * (y0 + (v0 + beta*y0)*t)

    else:
        # over-damped
        r1 = -beta + np.sqrt(beta**2 - omega0**2)
        r2 = -beta - np.sqrt(beta**2 - omega0**2)
        C1 = (v0 - r2*y0)/(r1 - r2)
        C2 = y0 - C1
        return C1*np.exp(r1*t) + C2*np.exp(r2*t)

plt.figure(figsize=(10, 6))

for i, values in enumerate(k):

    times = np.array(results['times'][values])
    y_num = np.array(results['y'][values])

    omega0 = np.sqrt(values/m)
    beta = gamma[i] / (2*m)

    y_exact = analytic_solution(times, y0, vy0, beta, omega0)
    error = np.abs(y_num - y_exact)

    plt.plot(times, error, color=colors[i],
             label=rf"$k={values}$ N/m")

plt.xlabel("Time $t$ (s)")
plt.ylabel(" $|y_\\mathrm{{num}} - y_\\mathrm{{exact}}|$")
plt.title("Numerical Error Rk4")
plt.yscale("log")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
#plt.savefig("\\numerical_error_spring.png", dpi=300, bbox_inches='tight')
plt.show()