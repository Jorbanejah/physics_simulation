import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# Para que el movimiento de la trayectoria se ve suave cuando cambie de medio,
# tenemos que hacer un metodo de incremental (pasando a la formulas diferenciales) con un tiempo t = 0.1 que sentara las bases del siguiente para x,y

class AnimatedProjectileMotion:

    g = 9.81
    k_values = np.linspace(0, 2, 20) 

    def __init__(self, x0, grades, v0):

        self.x = x0
        self.y = x0

        alpha = grades * np.pi / 180
        self.vx = v0 * np.cos(alpha)
        self.vy = v0 * np.sin(alpha)

        self.dt = 0.1
        self.frame = 0
        self.current_medium = 0 # Resistant medium
        self.boundaries = [] 

        self.fig, self.ax = plt.subplots() #return a tuple where: self.fig stores the figure object and self.ax stores the axes object
        self.point, = self.ax.plot([], [], marker='o', markersize=12, color='red')

        self.ax.set_xlim(0, 1000)
        self.ax.set_ylim(0, 300)
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")
        self.ax.set_title("Proyectil atravesando distintos medios")
        
    def get_medium(self,frame):
        return frame//20
    
    def update(self, frame):

        #New medium?
        new_medium = self.get_medium(frame)

        if new_medium != self.current_medium: 
            line = self.ax.axvline(self.x, color='gray', linestyle='--',) 
            self.boundaries.append(line)
            self.ax.text(
                self.x,               
                self.ax.get_ylim()[1] * frame/100,   # 95% height
                f"k = {self.k_values[new_medium]:.2f}",
                ha='center',
                va='top',
                fontsize=10,
                color='gray'
            )

            self.current_medium = new_medium
    
        k = self.k_values[new_medium]

        # Current velocity
        self.vx += - k * self.vx * self.dt
        self.vy += - self.g * self.dt - k * self.vy * self.dt

        # Current position
        self.x += self.vx * self.dt
        self.y += self.vy * self.dt

        # y >= 0
        if self.y < 0:
            self.y = 0
            raise StopIteration

        self.point.set_data([self.x], [self.y])
        return self.point,

    def animate(self):

        self.anim = FuncAnimation(self.fig, self.update, frames=100, interval=50, blit=False)
        plt.show()

if __name__ == "__main__":
    projectile_motion = AnimatedProjectileMotion(x0=0, grades=25, v0=100)
    projectile_motion.animate()