from SHM import Spring
import numpy as np
import matplotlib.pyplot as plt

def draw(t, t_max, step_spring):
    Energies = []
    times = []
    y_axes = []
    vy_phase = []
  
    while t < t_max:
        y, vy = step_spring.compute_motion(t)

        # Store the results for searching for crossings later
        y_axes.append(y)
        times.append(t)
        vy_phase.append(vy)
        #Energies
        _, _, Em = step_spring.compute_energy(y, vy)
        Energies.append(Em)

        t += step_spring.dt

    crossings_rk4 = get_period(y_axes, times, step_spring)
    
    return max(Energies), min(Energies), y_axes, vy_phase, crossings_rk4

def get_period(y, times, step_spring):
    crossings = []
    for i in range(1, len(y)):
        if y[i-1] < step_spring.y0 and y[i] >= step_spring.y0:  # ascending crossing
            #Linear interpolation to find accurate crossing time

            t1, y1 = times[i-1], y[i-1]
            t2, y2 = times[i], y[i]
            frac = -y1 / (y2 - y1)
            t_cross = t1 + frac * (t2 - t1)
            crossings.append(t_cross)

    if len(crossings) < 2:
        return None  # Not enough crossings to determine a period
    return crossings[1] - crossings[0] #Period

def T_teorico(k, m):
    T0 = 2 * np.pi * np.sqrt(m/k)
    return T0  

k = np.linspace(1, 20, 20)
Energy_average = np.zeros((2, len(k)))
T = []
T_teo = []
time = np.arange(0, 20, 0.05)

phase_space = {
    'k': k,
    'y': {},
    'vy': {}
}


for i, ang in enumerate(k):
    t_max = 20
    step_spring = Spring(k = ang, m = 1, A = 1, delta = 0, y0 = -2, x0 = 0, Anim = False, t_max = 20, dt=0.05)
    
    max_E, min_E, y_axes, vy_phase, crossings = draw(0, t_max, step_spring)

    T.append(crossings)
    T_teo.append(T_teorico(ang, step_spring.m))

    phase_space['y'][ang] = y_axes
    phase_space['vy'][ang] = vy_phase

    Energy_average[0][i] = max_E
    Energy_average[1][i] = min_E

fig1 = plt.figure(figsize=(10, 5))
plt.plot(k, Energy_average[0], label="Max Energy")
plt.plot(k, Energy_average[1], label="Min Energy")
plt.xlabel("k (N/m)")
plt.ylabel("Energy (J)")
plt.title("Energy vs k")

fig2 = plt.figure(figsize=(10, 5))
plt.plot(time, phase_space['y'][1.0], label = 'k = 1')
plt.plot(time, phase_space['y'][2.0], label = 'k = 2')
plt.xlabel('Time (s)')
plt.ylabel('y (m)')
plt.title('Spring motion for different values of k')
plt.legend()

fig3 = plt.figure(figsize=(10, 5))
plt.plot(k, T, label="Period")
plt.plot(k, T_teo, label="Theorical period")
plt.xlabel("k (N/m)")
plt.ylabel("Period (s)")
plt.title("Periods vs k")
plt.legend()    

fig4, ax = plt.subplots(2, 2, figsize=(10, 5))

j = 0
for i, ang in enumerate(k):
    if i % 5 == 0 and j < 4:
        ax.flat[j].plot(phase_space['y'][ang], phase_space['vy'][ang])
        ax.flat[j].set_xlabel("Velocity (m/s)")
        ax.flat[j].set_ylabel("y (m)")
        ax.flat[j].set_title(f"Phase space (k = {ang:.1f} N/m)")
        ax.flat[j].axhline(0, color='black', lw=0.5)
        ax.flat[j].axvline(step_spring.y0, color='black', lw=0.5)
        ax.flat[j].grid(True, alpha=0.3)
        j += 1

fig4.tight_layout()


cmap = plt.get_cmap('viridis')

fig5 = plt.figure(figsize=(8, 6))

# Normalize the spring constant k
k_norm = (k - np.min(k)) / (np.max(k) - np.min(k))

for ang, norm in zip(k, k_norm):
    color = cmap(norm)
    plt.plot(
        phase_space['y'][ang],
        phase_space['vy'][ang],
        color=color,
        linewidth=1
    )

sm = plt.cm.ScalarMappable(cmap=cmap, 
                           norm=plt.Normalize(vmin=np.min(k),
                                              vmax=np.max(k)))
sm.set_array([])

cbar = plt.colorbar(sm, ax=plt.gca()) 
cbar.set_label("k (N/m)")

plt.xlabel("y (m)")
plt.ylabel("vy (m/s)")
plt.title("Phase space with continuous colormap")
plt.grid(alpha=0.3)

plt.show()