import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

g = 9.81
x0, y0, grades, v0 = 0, 0, 45, 30 # Initial condition
t = np.linspace(0, 10, 100) #Time
alpha = grades*np.pi/ 180
k = [0, 0.05, 0.1, 0.2, 0.5]
v0x = v0*np.cos(alpha)
v0y = v0*np.sin(alpha)
x = np.zeros((len(t), len(k)))
y = np.zeros((len(t), len(k)))

#In this kind of motion the acceleration parameterer is function of velocity through this equations:
for i in range(len(k)):
    for j in range(len(t)):
        if k[i] == 0: # parabolic motion
            x[j][i] = v0x * t[j] + x0
            y[j][i] = -0.5 * g*t[j]**2 + v0y*t[j] + y0
        else:
            x[j][i] = v0x/k[i]*(1 - np.exp( - k[i]*t[j])) + x0
            y[j][i] = -g*t[j]/k[i] + (k[i]*v0y+g)/k[i]**2 * (1 - np.exp(-k[i]*t[j])) + y0

plt.figure()

for i in range(len(k)):
    plt.plot(x[:,i], y[:,i], label=f'k={k[i]}')

plt.grid()
plt.xlim(0, max(x.max(), 0))
plt.ylim(0, max(y.max(), 50))
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Projectile motion along several coefficient')
plt.legend()
plt.show()

#                                                    Numerical error (Euler's method) - animation.py
dt = 0.1

x_approx = [x0]
y_approx = [y0]
vx = v0x
vy = v0y

for i in range(len(k)):
    for j in range(len(t)):
        vx += - k[i] * vx * dt
        vy += - g * dt - k[i] * vy * dt

        # Current position
        x_approx.append(x_approx[-1] + vx * dt)
        y_approx.append(y_approx[-1] + vy * dt)



#                                            Analitical aproximations for range and time of flight
#Aprox k --> min

T_aprox = []
R_aprox = []
T = []
R = []
k1 = list(np.linspace(0, 0.1, 10))
for k in k1:
    if k == 0:
        T_aprox.append(2*v0y/g)
        R_aprox.append(v0x*(2*v0y/g))
        T.append(2*v0y/g)
        R.append(v0x*(2*v0y/g))
    else:
        T_aprox.append(2*v0y/g*(1 - k*v0y / (3*g)))
        R_aprox.append(v0x*(T_aprox[-1] - 0.5*k*T_aprox[-1]**2))

        #Solving transcendental equation for T numerically
        func = lambda t: t - (k*v0y+g)/(g*k)*(1 - np.exp(-k*t))
        T_real = fsolve(func, T_aprox[-1])[0]
        T.append(T_real)
        R.append(v0x*T_real)

plt.figure()
plt.plot(k1, T, label='Real time')
plt.plot(k1, T_aprox, label='Approximated time')
plt.xlabel('k')
plt.ylabel('Time (s)')
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(k1, R, label='Real range')
plt.plot(k1, R_aprox, label='Approximated range')
plt.xlabel('k')
plt.ylabel('Range (m)')
plt.legend()
plt.grid()
plt.show()