from SHM import Pendulum
import numpy as np
import matplotlib.pyplot as plt

def draw(t, t_max, pend_rk4, pend_apx):
    errors = []
    Energies_rk4 = []
    Energy_apx = []
    times = []
    thetas_rk4 = []
  

    while t < t_max:
        theta1, omega1 = pend_rk4.step_rk4()
        theta2, omega2 = pend_apx.step_approx()

        # Store the results for searching for crossings later
        thetas_rk4.append(theta1)
        times.append(t)

        #Energies
        _, _, Em_rk = pend_rk4.compute_energy(theta1, omega1)
        Energies_rk4.append(Em_rk)

        _, _, Em_aprx = pend_apx.compute_energy(theta2, omega2)
        Energy_apx.append(Em_aprx)

        errors.append(abs(theta1 - theta2))

        t += pend_rk4.dt

    crossings_rk4 = get_period(thetas_rk4, times)
    
    return np.mean(errors), np.std(errors), max(Energies_rk4), min(Energies_rk4), max(Energy_apx), min(Energy_apx), crossings_rk4

def get_period(thetas, times):
    crossings = []
    for i in range(1, len(thetas)):
        if thetas[i-1] < 0 and thetas[i] >= 0:  # ascending crossing
            #Linear interpolation to find accurate crossing time

            t1, th1 = times[i-1], thetas[i-1]
            t2, th2 = times[i], thetas[i]
            frac = -th1 / (th2 - th1)
            t_cross = t1 + frac * (t2 - t1)
            crossings.append(t_cross)

    if len(crossings) < 2:
        return None  # Not enough crossings to determine a period
    return crossings[1] - crossings[0]

def T_teorico(theta0, L=1, g=9.81):
    T0 = 2 * np.pi * np.sqrt(L / g)
    return T0 * (1 + (1/16)*theta0**2 + (11/3072)*theta0**4) # Taylor expansion

errors = []
theta = np.linspace(0, np.deg2rad(90), 50)
Energy_average = np.zeros((4, len(theta)))
t_max = 20
T_rk4 = []
T_teo = []

for i, ang in enumerate(theta):
    pend_rk4 = Pendulum(L=1, theta0=ang, omega0=0, m=1, t_max=t_max, dt=0.01)
    pend_apx = Pendulum(L=1, theta0=ang, omega0=0, m=1, t_max=t_max, dt=0.01)

    error, sigma, max_E_rk4, min_E_rk4, max_E_apx, min_E_apx, crossings_rk4 = draw(0, t_max, pend_rk4, pend_apx)

    T_rk4.append(crossings_rk4)
    T_teo.append(T_teorico(ang))
   # T_apx.append(get_period(pend_apx, t_max))

    Energy_average[0][i] = max_E_rk4
    Energy_average[1][i] = min_E_rk4
    Energy_average[2][i] = max_E_apx
    Energy_average[3][i] = min_E_apx

    errors.append(error)

# Until where the approximation is valid, the energy should be constant and equal to the initial energy.
# We can see how the energy starts to deviate as the initial angle increases, showing that the approximation is no longer valid. 
fig1 = plt.figure(figsize=(10, 5))
plt.plot(np.rad2deg(theta), Energy_average[2], label="Max Energy Apx")
plt.plot(np.rad2deg(theta), Energy_average[3], label="Min Energy Apx")
plt.xlabel("Initial Angle (degrees)")
plt.ylabel("Energy (J)")
plt.title("Energy vs Initial Angle")
plt.legend()

fig2 = plt.figure(figsize=(10, 5))
plt.plot(np.rad2deg(theta), Energy_average[0], label="Max Energy RK4")
plt.plot(np.rad2deg(theta), Energy_average[1], label="Min Energy RK4")
plt.xlabel("Initial Angle (degrees)")
plt.ylabel("Energy (J)")
plt.title("Energy vs Initial Angle")
plt.legend()

# The error should be small for small angles and increase as the angle increases, showing that the apoorximation is no longer valid.
fig3 = plt.figure(figsize=(10, 5))
plt.plot(np.rad2deg(theta), errors)
plt.xlabel('Initial Angle (degrees)')
plt.ylabel('Energy Error')
plt.title('Error vs Initial Angle')

fig4 = plt.figure(figsize=(10, 5))
plt.plot(np.rad2deg(theta), T_rk4, label="Periodo RK4")
plt.plot(np.rad2deg(theta), T_teo, label="Periodo Teórico")
plt.xlabel("Initial Angle (degrees)")
plt.ylabel("Period (s)")
plt.title("Period vs Initial Angle")
plt.legend()    

plt.show()

