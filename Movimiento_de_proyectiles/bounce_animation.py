import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Bounce_animation:

    g = 9.8

    def __init__(self, x0, y0, v0, angle, k, e):

        self.x = x0
        self.y = y0
        self.v0 = v0
        alpha = np.radians(angle)
        self.vx = v0 * np.cos(alpha)
        self.vy = v0 * np.sin(alpha)

        self.k = k #Resistance parameter
        self.e = e #Bounce coefficient Restitution parameter

        self.dt = 0.1 # Time step

        self.fig, self.ax = plt.subplots()
        self.point, = self.ax.plot([], [], marker='o', markersize=12, color='blue')

        H = self.v0 **2 * np.sin(alpha)**2 / (2*self.g)
        R = self.v0 **2 * np.sin(2*alpha) / self.g
        if self.e == 0:
            self.ax.set_xlim(0, R + 10)
            self.ax.set_ylim(0, H + 10)
        else:
            self.ax.set_xlim(0, R + R/2)
            self.ax.set_ylim(0, H + 10)

        #Graphics features
        self.ax.set_xlabel('x (m)')
        self.ax.set_ylabel('y (m)')
        self.ax.set_title(f"Bounce coefficient e = {self.e}" f"Resistance parameter k = {self.k}")

    def update(self, frame):
        
        #current velocity
        self.vx += - self.k * self.vx * self.dt
        self.vy += - self.g * self.dt - self.k * self.vy * self.dt

        # Current position
        self.x += self.vx * self.dt
        self.y += self.vy * self.dt

        if self.y <= 0:
            self.y = 0
            self.vy = - self.vy * self.e

            if self.vy < 1e-1:
                self.anim.event_source.stop()

            
        self.point.set_data([self.x + 6], [self.y + 6])
        return self.point,

    def animate(self):
        self.animate = FuncAnimation(self.fig, self.update, frames =100, interval = 50, blit = False)
        
        plt.show()

if __name__ == '__main__':
    projectile_motion = Bounce_animation(0, 0, 100, 45, 0.2, 1)
    projectile_motion.animate()