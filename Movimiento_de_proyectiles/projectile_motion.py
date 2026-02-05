import numpy as np
import matplotlib.pyplot as plt

#################################################   Initial parameters   ###################

g = 9.81 #Gravity
x0, y0, grades, v0 = 0, 0, 45, 30 # Initial condition
t = np.linspace(0, 10, 100) #Time
alpha = np.radians(grades)
dt = t[1] - t[0] #Step
k_value = [0, 0.05, 0.1, 0.2, 0.5]
v0x = v0*np.cos(alpha)
v0y = v0*np.sin(alpha)
e = 0.8 #Rebote 0 < e < 1

################################################# Movement of projectiles through different media

#  Analitical
x = np.zeros((len(t), len(k_value)))
y = np.zeros((len(t), len(k_value)))
vx = []
vy = []

for i, k in enumerate(k_value): 
    if k == 0: 
        x[:, i] = v0x * t + x0 
        y[:, i] = -0.5*g*t**2 + v0y*t + y0 
        vx.append(v0x)
        vy.append(v0y - g*t)
        y_col = y[:, i]
        idx_ground = np.where(y_col < 0)[0]

        if len(idx_ground) > 0:
            first_hit = idx_ground[0]
            y_col[first_hit:] = 0
            continue

    else: 
        x[:, i] = v0x/k * (1 - np.exp(-k*t)) + x0 
        y[:, i] = -g*t/k + (k*v0y + g)/k**2 * (1 - np.exp(-k*t)) + y0
        vx.append(v0x * np.exp( - k * t))
        vy.append(-g/k + (k*v0y + g)/k * np.exp(- k * t))

        y_col = y[:, i]
        idx_ground = np.where(y_col < 0)[0]

        if len(idx_ground) > 0:
            first_hit = idx_ground[0]
            y_col[first_hit:] = 0
            continue

# Numerical Euler's method

x_numeric = np.zeros((len(t), len(k_value)))
y_numeric = np.zeros((len(t), len(k_value)))
vx_numeric = []
vy_numeric = []

for i, k in enumerate(k_value):
    x_pos = x0
    y_pos = y0
    vx = v0x
    vy = v0y
    for j in range(len(t)):

        x_numeric[j][i] = x_pos
        y_numeric[j][i] = y_pos
        vx_numeric.append(vx)
        vy_numeric.append(vy)

        #Current velocity
        vx += - k * vx * dt
        vy += - g * dt - k * vy * dt
        # Current position
        x_pos += vx * dt
        y_pos += vy * dt

        if y_pos < 0:
            y_pos = 0
            continue

#  Plot

fig1, axs1 = plt.subplots(2, 2, figsize = (10, 8))
axs1 = axs1.flatten()
for idx, k in enumerate(k_value[1:]):   
    col = idx + 1 

    axs1[idx].plot(x[:,col], y[:,col], label="Exacta")
    axs1[idx].plot(x_numeric[:, col], y_numeric[:, col], '--', label="Euler")

    axs1[idx].set_title(f"k = {k}")
    axs1[idx].set_xlabel("x (m)")
    axs1[idx].set_ylabel("y (m)")
    axs1[idx].grid()
    axs1[idx].legend()

plt.tight_layout()

################################################ Rebound movement

# Analitical
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
  

# Numerical

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