import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

g = 9.81
x0, y0, grades, v0 = 0, 0, 45, 30 # Initial condition
t = np.linspace(0, 10, 100) #Time
alpha = grades*np.pi/ 180
k_value = [0, 0.05, 0.1, 0.2, 0.5]
v0x = v0*np.cos(alpha)
v0y = v0*np.sin(alpha)
x = np.zeros((len(t), len(k_value)))
y = np.zeros((len(t), len(k_value)))

e = 0.8 #Rebote 0 < e < 1

#In this kind of motion the acceleration parameterer is function of velocity through this equations:


# al final de cada k
x_totales =[]
y_totales =[]

for i, k in enumerate(k_value):
    x_posi = x0
    y_posi = y0
    vx_posi = v0x
    vy_posi = v0y

    t_local = t.copy()

    x_list_rebot = []
    y_list_rebot = []
    maximos_rebotes = 100

    for j in range(maximos_rebotes):

        if k == 0:
            x_tramo= vx_posi * t_local + x_posi 
            y_tramo = -0.5*g*t_local**2 + vy_posi*t_local + y_posi

        else:
            x_tramo= vx_posi/k * (1 - np.exp(-k*t_local)) + x_posi 
            y_tramo= -g*t_local/k + (k*vy_posi + g)/k**2 * (1 - np.exp(-k*t_local)) + y_posi

        idx = np.where(y_tramo < 0 )[0]

        if len(idx) == 0:
            x_list_rebot.append(x_tramo)
            y_list_rebot.append(y_tramo)
            break
        
        impact = idx[0]

        x_impact = x_tramo[impact]
        t_impact = t_local[impact]

        if k == 0:
            vy_impact = vy_posi - g * t_impact
            vx_impact = vx_posi
        else:
            vx_impact = vx_posi * np.exp( - k * t_impact)
            vy_impact = -g/k + (k*vy_posi + g)/k * np.exp(- k * t_impact)


        x_list_rebot.append(x_tramo[:impact + 1]) 
        y_list_rebot.append(y_tramo[:impact + 1])

        #Update current position

        x_posi = x_impact
        y_posi = 0
        vy_posi = - vy_impact * e
        vx_posi = vx_impact
        
        if abs(vy_posi) < 1e-2 and abs(vx_posi) < 1e-2:
            break

        t_local = t_local - t_impact
        t_local = t_local[t_local >= 0]

    
    x_totales.append(np.concatenate(x_list_rebot))
    y_totales.append(np.concatenate(y_list_rebot))



        
for i, k in enumerate(k_value): 
    if k == 0: 
        x[:, i] = v0x * t + x0 
        y[:, i] = -0.5*g*t**2 + v0y*t + y0 

        y_col = y[:, i]
        idx_ground = np.where(y_col < 0)[0]

        if len(idx_ground) > 0:
            first_hit = idx_ground[0]
            y_col[first_hit:] = 0

    else: 
        x[:, i] = v0x/k * (1 - np.exp(-k*t)) + x0 
        y[:, i] = -g*t/k + (k*v0y + g)/k**2 * (1 - np.exp(-k*t)) + y0#

        y_col = y[:, i]
        idx_ground = np.where(y_col < 0)[0]

        if len(idx_ground) > 0:
            first_hit = idx_ground[0]
            y_col[first_hit:] = 0
#                                                    Numerical error (Euler's method) - ref animation.py

dt = t[1] - t[0]

x_approx = np.zeros((len(t), len(k_value)))
y_approx = np.zeros((len(t), len(k_value)))

for i, k in enumerate(k_value):
    x_pos = x0
    y_pos = y0
    vx = v0x
    vy = v0y
    for j in range(len(t)):

        x_approx[j][i] = x_pos
        y_approx[j][i] = y_pos

        #Current velocity
        vx += - k * vx * dt
        vy += - g * dt - k * vy * dt
        # Current position
        x_pos += vx * dt
        y_pos += vy * dt

        if y_pos < 0:
            y_pos = 0
            vy = -e * vy
        if y_pos == 0 and abs(vy)<1e-2:
            vy = 0
            break


fig, axs = plt.subplots(2,2, figsize=(10, 8))
axs = axs.flatten()

for idx, k in enumerate(k_value[1:]):   
    col = idx + 1 
    x_total = x_totales[col]
    y_total = y_totales[col]

    mask_exact = y_total >= 0
    mask_num   = y_approx[:, col] >= 0

    axs[idx].plot(x_total[mask_exact], y_total[mask_exact], label="Exacta")
    axs[idx].plot(x_approx[mask_num, col], y_approx[mask_num, col], '--', label="Euler")

    axs[idx].set_title(f"k = {k}")
    axs[idx].set_xlabel("x (m)")
    axs[idx].set_ylabel("y (m)")
    axs[idx].grid()
    axs[idx].legend()

plt.tight_layout()
plt.show()


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