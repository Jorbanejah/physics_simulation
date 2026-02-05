import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve

g = 9.81 #Gravity
x0, y0, grades, v0 = 0, 0, 45, 30 # Initial condition
alpha = np.radians(grades)
v0x = v0*np.cos(alpha)
v0y = v0*np.sin(alpha)
m = 0.8

def energy(m, g, x, y, vx, vy, E0=None):
    """
    Calcula energías para cada k:
    - E_mec: total mechanic energy
    - E_k: kinetic energy
    - E_p: potencial energy
    - P_ext: losses
    """
    x = np.array(x)
    y = np.array(y)
    vx = np.array(vx)
    vy = np.array(vy)

    Ek = 0.5 * m * (vx**2 + vy**2)
    Ep = m * g * y 
    Em = Ek + Ep

    if E0 is None:
        E0 = Em[0] if hasattr(Em, '__len__') else Em

    P_ext = E0 - Em
    return Em, Ek, Ep, P_ext

# Time array
t = np.linspace(0, 10, 100) #Time
dt = t[1] - t[0] #Step
k_value = [0, 0.1]

# Calculate analytical and numerical solutions for each k
results = {}

for k in k_value:
    # Analytical method

    x_ana = v0x/k * (1 - np.exp(-k*t)) + x0 if k != 0 else v0x * t + x0
    y_ana = (-g*t/k + (k*v0y + g)/k**2 * (1 - np.exp(-k*t)) + y0) if k != 0 else (-0.5*g*t**2 + v0y*t + y0)
    vx_ana = v0x * np.exp(-k*t) if k != 0 else np.full_like(t, v0x)
    vy_ana = (-g/k + (k*v0y + g)/k * np.exp(-k*t)) if k != 0 else (v0y - g*t)
    
    # Numerical method (Euler)
    x_num, y_num, vx_num, vy_num = [x0], [y0], [v0x], [v0y]
    
    for _ in range(len(t)-1):
        vx_new = vx_num[-1] - k * vx_num[-1] * dt
        vy_new = vy_num[-1] - g * dt - k * vy_num[-1] * dt
        x_new = x_num[-1] + vx_new * dt
        y_new = y_num[-1] + vy_new * dt
        
        if y_new < 0:
            y_new = 0
        
        x_num.append(x_new)
        y_num.append(y_new)
        vx_num.append(vx_new)
        vy_num.append(vy_new)
    
    # Calculate energies
    Em_ana, Ek_ana, Ep_ana, P_ana = energy(m, g, x_ana, y_ana, vx_ana, vy_ana)
    Em_num, Ek_num, Ep_num, P_num = energy(m, g, x_num, y_num, vx_num, vy_num)
    
    results[k] = {
        'Em_ana': Em_ana, 'Em_num': Em_num,
        'Ek_ana': Ek_ana, 'Ek_num': Ek_num,
        'Ep_ana': Ep_ana, 'Ep_num': Ep_num,
        'P_ana': P_ana, 'P_num': P_num
    }
    # Plot energies over time
    
fig, axes = plt.subplots(1, 1, figsize=(14, 5))
axes = axes.flatten()
for idx, k in enumerate(k_value):
    axes[idx].plot(t, results[k]['Ek_ana'], label='Ek (analytical)', linewidth=2)
    axes[idx].plot(t, results[k]['Ek_num'], '--', label='Ek (numerical)', linewidth=2)
    axes[idx].plot(t, results[k]['Ep_ana'], label='Ep (analytical)', linewidth=2)
    axes[idx].plot(t, results[k]['Ep_num'], '--', label='Ep (numerical)', linewidth=2)    

    axes[idx].set_xlabel('Time (s)')
    axes[idx].set_ylabel('Energy (J)')
    axes[idx].set_title(f'Energy vs Time (k = {k})')
    axes[idx].legend()
    axes[idx].axhline(0, color='gray', lw=0.5, ls='--')


plt.tight_layout()
    
fig1, ax1 = plt.subplots(1, 1, figsize=(7, 5))

for idx, k in enumerate(k_value):
    ax1.plot(t, results[k]['Em_ana'], label='Em (analytical)', linewidth=2)
    ax1.plot(t, results[k]['Em_num'], '--', label='Em (numerical)', linewidth=2)
    ax1.plot(t, results[k]['P_ana'], label='Losses', linewidth=2)
    ax1.plot(t, results[k]['P_num'], '--', label='Losses (numerical)', linewidth=2)

    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Energy (J)')
    ax1.set_title(f'Mechanical Energy and Losses vs Time (k = {k})')
    ax1.legend()
    ax1.axhline(0, color='gray', lw=0.5, ls='--')

plt.tight_layout()

#                                                 Analitical method - Range, Height and time
# Everybody knows how to get the range and height of a parabolic motion. Through the parabolic equation when y = 0 we can solve it and see that 
# R = x = self.v0 **2 * np.sin(2*alpha) / self.g and then we can get the time. On the other hand, when velocity vy = 0 we can solve the equation put inside y equation and take it, too: H = self.v0 **2 * np.sin(alpha)**2 / (2*self.g).
# However, now we are trying to figure out how to deal with parabolic equation through different medium. So, putting y = 0 we can see a trascental equation that we can solve by numerical method (for example 
# plotting both parts right and left, and show where interline) or doing some approximation such as exponential Taylor expansion.


T_aprox = [] # Approx time
R_aprox = [] # Approx range


T = []
R = []
H = []
k1 = list(np.linspace(0, 0.8, 20))

for k in k1:
    
    if k == 0:
        T_height = v0y/g
        T_aprox.append(2*v0y/g)
        R_aprox.append(v0x*(2*v0y/g))
        
        H.append(y0 + v0y*T_height - 1/2*g*T_height**2)
        T.append(2*v0y/g)
        R.append(v0x*(2*v0y/g))

    else:
        T_height = -1/k*np.log(g/(k*v0y+g))

        T_aprox.append(2*v0y/g*(1 - k*v0y / (3*g)))
        R_aprox.append(v0x*(T_aprox[-1] - 0.5*k*T_aprox[-1]**2))

        #Solving transcendental equation for T -- numerically
        func = lambda t: t - (k*v0y+g)/(g*k)*(1 - np.exp(-k*t))
        T_real = fsolve(func, T_aprox[-1])[0]
        T.append(T_real)
        R.append(v0x/k*(1- np.exp(-k*T_real)))
        H.append((k*v0y + g)/k**2 * (1 - np.exp(-k*T_height)) - g*T_height/k)



plt.figure(1)
plt.plot(k1, T, label='Real time')
plt.plot(k1, T_aprox, '--',label='Approximated time')
plt.title('Time in air vs resistance')
plt.xlabel('k (1/s)')
plt.ylabel('Time (s)')
plt.legend()
plt.show()


plt.figure(2)
plt.plot(k1, R, label='Real range')
plt.plot(k1, R_aprox, '--',label='Approximated range')
plt.title('Maximum range vs resistance')
plt.xlabel('k (1/s)')
plt.ylabel('Range (m)')
plt.legend()
plt.show()


plt.figure(3)
plt.plot(k1, H, '*',label='Height')
plt.xlabel('k (1/s)')
plt.ylabel('Height (m)')
plt.title('Maximum height vs resistance')
plt.legend()
plt.show()