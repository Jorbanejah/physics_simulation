
import numpy as np
from Damping_vibration import Spring
import matplotlib.pyplot as plt


def draw(t, t_max, spr):
    Ep = []
    Ek = []
    Wp = []
    times = []
    y_axes = []
    vy_axes= []
  
    while t < t_max:
        y, vy = spr.rk4(t)

        # Store the results for searching for crossings later
        y_axes.append(y)
        times.append(t)
        vy_axes.append(vy)

        #Energies
        Ek0, Ep0, Wp0, _ = spr.energy(y, vy, t)
        Ep.append(Ep0)
        Ek.append(Ek0)
        Wp.append(Wp0)

        t += spr.dt

    Et = np.array(Ep) + np.array(Ek) + np.array(Wp)

    return Ep, Ek, Wp, Et, times, y_axes, vy_axes


k = [2, 8, 18]
m = 2
gamma = [2, 8, 13]

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
    spr = Spring(-2, 0, m, values, gamma[i], dt = 0.01, t_max =15)

    Ep, Ek, Wp, Et, times, y, vy = draw(t, t_max, spr)

    results['y'][values] = y
    results['vy'][values] = vy
    results['Ep'][values] = Ep
    results['Ek'][values] = Ek
    results['Wp'][values] = Wp
    results['Et'][values] = Et
    results['times'][values] = times

##
# --------- Graphics ---------
##

# y vs t with different k

figure1 = plt.figure(figsize=(10, 6))
for i, values in enumerate(k):
    plt.plot(results['times'][values], results['y'][values], label = f'k={values:.1f} (N/m)')
plt.xlabel('Time (t)')
plt.ylabel('y (m)')
plt.legend()
plt.show()

figure2, ax2 = plt.subplots(1, 3, figsize= (14,8), tight_layout = True)
for i, values in enumerate(k):
    ax2[i].plot(results['times'][values], results['Ek'][values], label = 'Ek')
    ax2[i].plot(results['times'][values], results['Ep'][values], label = 'Ep')
    ax2[i].plot(results['times'][values], results['Wp'][values], label = 'Wp')
    ax2[i].plot(results['times'][values], results['Et'][values], label = 'Et')

    ax2[i].set_xlabel('Time (s)')
    ax2[i].set_ylabel('Energy (J)')
    ax2[i].set_title(f'Energy k = {values:.1f} (N/m)')
    ax2[i].legend()


plt.show()
##
# --------- Colormaps ---------
##
##
# --------- Errors RK4 --------
##