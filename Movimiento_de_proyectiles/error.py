import numpy as np
import matplotlib.pyplot as plt

#                                                  Analitical aproximations for range and time of flight
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

        #Solving transcendental equation for T -- numerically
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