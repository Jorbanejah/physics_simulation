from SHM import Spring
import numpy as np
import matplotlib.pyplot as plt

def draw(t, t_max, step_spring):
    Energies = []
    times = []
    y_axes = []
  
    while t < t_max:
        y, vy = step_spring.compute_motion(t)

        # Store the results for searching for crossings later
        y_axes.append(y)
        times.append(t)

        #Energies
        _, _, Em = step_spring.compute_energy(y, vy)
        Energies.append(Em)

        t += step_spring.dt

    crossings_rk4 = get_period(y_axes, times, step_spring)
    
    return max(Energies), min(Energies), y_axes, crossings_rk4

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
y_compute = np.zeros((len(k), len(time)))


for i, ang in enumerate(k):
    t_max = 20
    step_spring = Spring(k = ang, m = 1, A = 1, delta = 0, y0 = -2, x0 = 0, Anim = False, t_max = 20, dt=0.05)
    
    max_E, min_E, y_axes, crossings = draw(0, t_max, step_spring)

    T.append(crossings)
    T_teo.append(T_teorico(ang, step_spring.m))

    y_compute[i][:] = y_axes
    Energy_average[0][i] = max_E
    Energy_average[1][i] = min_E

fig1 = plt.figure(figsize=(10, 5))
plt.plot(k, Energy_average[0], label="Max Energy")
plt.plot(k, Energy_average[1], label="Min Energy")
plt.xlabel("k (N/m)")
plt.ylabel("Energy (J)")
plt.title("Energy vs k")


fig2 = plt.figure(figsize=(10, 5))
plt.plot(time, y_compute[0], label = 'k = 1')
plt.plot(time, y_compute[2], label = 'k = 2')
plt.legend(k)
plt.xlabel('Time (s)')
plt.ylabel('y (m)')
plt.title('Spring motion for different values of k')
plt.legend(k)

fig3 = plt.figure(figsize=(10, 5))
plt.plot(k, T, label="Period")
plt.plot(k, T_teo, label="Theorical period")
plt.xlabel("k (N/m)")
plt.ylabel("Period (s)")
plt.title("Periods vs k")
plt.legend()    

plt.show()