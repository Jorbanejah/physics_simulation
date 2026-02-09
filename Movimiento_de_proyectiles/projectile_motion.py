import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# =========================================================
# Initial parameters
# =========================================================

g = 9.81
x0, y0 = 0.0, 0.0
v0 = 30.0
angle = 45
alpha = np.radians(angle)

v0x = v0 * np.cos(alpha)
v0y = v0 * np.sin(alpha)

k_values = [0, 0.05, 0.1, 0.2, 0.5]
e = 0.8                      # coeficiente de restitución
t = np.linspace(0, 10, 1000)
dt = t[1] - t[0]

# =========================================================
# Analitical solution (whithout rebound)
# =========================================================

def analytic_motion(k):
    if k == 0:
        x = v0x * t + x0
        y = -0.5 * g * t**2 + v0y * t + y0
        vx = np.full_like(t, v0x)
        vy = v0y - g * t
    else:
        x = v0x / k * (1 - np.exp(-k * t)) + x0
        y = -(g / k) * t + (k * v0y + g) / k**2 * (1 - np.exp(-k * t)) + y0
        vx = v0x * np.exp(-k * t)
        vy = -g / k + (k * v0y + g) / k * np.exp(-k * t)

    impact = np.where(y < 0)[0]
    if len(impact) > 0:
        y[impact[0]:] = 0
        x[impact[0]:] = x[impact[0]]
        vx[impact[0]:] = 0
        vy[impact[0]:] = 0

    return x, y, vx, vy

# =========================================================
# Numerical solution (whithout rebound)
# =========================================================

def euler_motion(k):
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    vx = np.zeros_like(t)
    vy = np.zeros_like(t)

    x[0], y[0] = x0, y0
    vx[0], vy[0] = v0x, v0y

    for i in range(1, len(t)):
        vx[i] = vx[i-1] - k * vx[i-1] * dt
        vy[i] = vy[i-1] - g * dt - k * vy[i-1] * dt

        x[i] = x[i-1] + vx[i] * dt
        y[i] = y[i-1] + vy[i] * dt

        if y[i] < 0:
            y[i] = 0
            vx[i] = 0
            vy[i] = 0
            continue

    return x, y, vx, vy

# =========================================================
# Analitical solution (with rebound)
# =========================================================

def analytic_rebound(k):
    x_all, y_all = [], []
    vx_all, vy_all = [], []

    x_pos, y_pos = x0, y0
    vx_pos, vy_pos = v0x, v0y
    t_local = t.copy()

    for _ in range(50):
        if k == 0:
            x = vx_pos * t_local + x_pos
            y = -0.5 * g * t_local**2 + vy_pos * t_local + y_pos

            vx = np.full_like(t_local, vx_pos)
            vy = vy_pos - g * t_local
        else:
            x = vx_pos / k * (1 - np.exp(-k * t_local)) + x_pos
            y = -(g / k) * t_local + (k * vy_pos + g) / k**2 * (1 - np.exp(-k * t_local)) + y_pos

            vx = vx_pos * np.exp(-k * t_local)
            vy = -g / k + (k * vy_pos + g) / k * np.exp(-k * t_local)

        hit = np.where(y < 0)[0]
        if len(hit) == 0:
            x_all.append(x)
            y_all.append(y)
            vx_all.append(vx)
            vy_all.append(vy)
            break

        i = hit[0]

        x_all.append(x[:i+1])
        y_all.append(y[:i+1])
        vx_all.append(vx[:i+1])
        vy_all.append(vy[:i+1])

        x_impact = x[i]
        t_impact = t_local[i]

        if k ==0:
            vy_impact = vy_pos - g * t_impact
        else:
            vy_impact = -g / k + (k * vy_pos + g) / k * np.exp(-k * t_impact)

        vy_pos = -e * vy_impact
        vx_pos = vx[i]
        x_pos = x[i]
        y_pos = 0.0

        if abs(vy_pos) < 1e-2:
            break
        
        t_hit = t_local[i]
        t_local = t_local - t_hit
        t_local = t_local[t_local >= 0]

    return (np.concatenate(x_all), np.concatenate(y_all), np.concatenate(vx_all), np.concatenate(vy_all))

# =========================================================
# Numerical solution (with rebound)
# =========================================================

def euler_rebound(k):
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    vx = np.zeros_like(t)
    vy = np.zeros_like(t)
    x[0], y[0] = x0, y0
    vx[0], vy[0] = v0x, v0y

    for i in range(1, len(t)):
        vx[i] = vx[i-1] - k * vx[i-1] * dt
        vy[i] = vy[i-1] - g * dt - k * vy[i-1] * dt

        x[i] = x[i-1] + vx[i] * dt
        y[i] = y[i-1] + vy[i] * dt

        if y[i] < 0:
            y[i] = 0
            vy[i] = -e * vy[i-1]

        if y[i] == 0 and abs(vy[i]) < 1e-2:
            vy[:i+1] = 0
            continue

    return x, y, vx, vy

# =========================================================
# Printing results
# =========================================================

fig = plt.figure(figsize=(10, 6))
x, y, _, _ = analytic_motion(k_values[0])
plt.plot(x, y, label="Analitical")
plt.xlim(0, np.max(x) + 10)
plt.ylim(0, np.max(y) + 10)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.savefig("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\Movimiento_de_proyectiles\\projectile_motion_analytic.png", dpi=300, bbox_inches='tight')


fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs = axs.flatten()

for i, k in enumerate(k_values[1:]):
    xa, ya, _, _ = analytic_motion(k)
    xn, yn, _, _ = euler_motion(k)

    axs[i].plot(xa, ya, label="Analitical")
    axs[i].plot(xn, yn, "--", label="Euler")
    axs[i].set_title(f"k = {k}")
    axs[i].set_xlabel("x (m)")
    axs[i].set_ylabel("y (m)")
    axs[i].grid()
    axs[i].legend()

plt.tight_layout()
plt.savefig("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\Movimiento_de_proyectiles\\projectile_motion_comparison.png", dpi=300, bbox_inches='tight')

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs = axs.flatten()

for i, k in enumerate(k_values[1:]):
    xa, ya, _, _ = analytic_rebound(k)
    xn, yn, _, _ = euler_rebound(k)

    axs[i].plot(xa, ya, label="Analitical")
    axs[i].plot(xn, yn, "--", label="Euler's method")
    axs[i].set_title(f"k = {k}")
    axs[i].grid()
    axs[i].legend()

plt.tight_layout()
plt.savefig("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\Movimiento_de_proyectiles\\projectile_motion_rebound_comparison.png", dpi=300, bbox_inches='tight')

# =========================================================
# Exporting data to .xlsx
# =========================================================

data = {}

for k in k_values:
    xa, ya, vxa, vya = analytic_motion(k)
    xn, yn, vxn, vyn = euler_motion(k)

    data[f"x_a_k={k}"] = xa
    data[f"y_a_k={k}"] = ya
    data[f"vx_a_k={k}"] = vxa
    data[f"vy_a_k={k}"] = vya
    data[f"x_n_k={k}"] = xn
    data[f"y_n_k={k}"] = yn
    data[f"vx_n_k={k}"] = vxn
    data[f"vy_n_k={k}"] = vyn

data["t"] = t

df = pd.DataFrame(data)

wb = Workbook()
ws = wb.active

for row in dataframe_to_rows(df, index=False, header=True):
    ws.append(row)

wb.save("C:\\Users\\JORGE\\Desktop\\Programas\\Python\\physics_simulation\\Movimiento_de_proyectiles\\projectile_motion_data.xlsx")