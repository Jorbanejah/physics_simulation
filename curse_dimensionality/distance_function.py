'''
The following program contains a class calls as the file.

This class shows us how high-dimensional space makes bigger when the dimension increase. 

Given the dimension (d), the numbers of point (N) and the norm type, the following class will calculate the distance (norm type) between two random 
d-dimensional points N times and plot it in a two-dimensional graphic (N vs avergare distance between points).

The main purpose of this class are both curiosity and show one performance of Machine Learning: dimensionaly reduction.
'''

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import time
from typing import Sequence
from matplotlib.figure import Figure

def auditor(func):
    def wrapper(*arg, **kwargs):

        if not arg:
            raise ValueError("The function requieres an argument")
        if not isinstance(*arg[0], int):
            raise TypeError("The first argument must be an int")
        
        print(f"Log... calling {func.__name__} with arg = {arg} and kwargs = {kwargs}")

        start = time.time()
        score = func(*arg, **kwargs)
        end = time.time()

        duracion = (start - end) *1000

        print(f"Score: {score}")

        return score
    return wrapper


class Hypercircle():
    """
    Description
    ------------

    This class calculates the percentage of p-dimensional data points that fall inside a hypersphere with a given center and radius. The data is generated from a normal distribution with a given mean and s.t.d.

    The class works as follows: first, the data is generated from p-dimensional normal distribution, and transformed into a dataset of size (n x p), called X. Finally, the proportion of points lying within the hypersphere is computed by a function called: in_or_out.

    Parameters:
    -----------

    p: int
        Dimension of the whole dataset
    n: int
        number of samples
    center: Array
        center of the hypercircle
    ratio: int
        ratio of the hypercircle 
    graphics_2d: bool
        the class generated a scatter picture of two dimensional data with a circle. 
    graphics_3d: bool
        the class generated a scatter picture of three dimensional data with a sphere.
    pd-graphcis: bool
        the class generate a curve about how the porcentage of points inside the hypercircle decrease with dimesion.
    """

    def __init__(self, p: int, n: int, center: list, radius: int, graphics_2d: bool = False, graphics_3d: bool= False, graphics_pd: bool = True, **kwargs):

        self.p: int = p
        self.n: int = n
        self.center: list = center
        self.radius:int = radius
        self.graphics_pd:bool = graphics_pd
        self.graphics_3d:bool = graphics_3d
        self.graphics_2d: bool = graphics_2d

        self.kwargs = kwargs

    def generate_data(self, p:int, n:int, ):
        """
        This function generate a DataFrame with the whole generate p-dimensional data given a mean and std. 
        By default, the mean and the std are 0 and 1, respectively.
        """

        mean = 0 if self.kwargs["mean"] is None else mean = self.kwargs["mean"]
        scale = 1 if self.kwargs["std"] is None else scale = self.kwargs["std"]

        try: 
            x0 = np.random.normal(loc = mean, scale = scale, size = (p, n))
            self.df = pd.DataFrame(x0.reshape(n, p))

            return self.df
        except:
            raise ValueError("The mean and the standard desviation (std) must be integer")

    def stratified(self, df, p:int) ->list[float]:

        self.porcentage = []

        for dimension in range(p + 1):

            if dimension == 0: continue

            new_df = df[:p][:]

            self.porcentage.append(self.porcentage_inside(df = new_df, p = dimension, n = self.n))

        return self.porcentage

    def porcentage_inside(self, df, p: int, n:int,)-> float:
        """
        This function calculates the porcentage of p-dimensional points that it falls inside the hypersphere with given radius and p-dimensional center.
        """
        j = 0

        for i in range(n):
            #We make a p-dimensional list for the whole points in the same column DataFrame

            point = [df[d][i] for d in range(p)]
            
            center_coordenate = [self.center[d] for d in range(p)]

            if self.in_or_out(self, center = center_coordenate, radio = self.radius, point = point) == True:
                j += 1

        porcent = j/n * 100

        return porcent

    def in_or_out(self, center:list, radio:int, point:list) -> bool:
        """
        This boolean function tells you whether a point is inside the hypersphere with a given radius and center or not.
        """
        center= np.array(center)
        point = np.array(point)
        distance = np.linalg.norm(center - point)

        if distance > radio:
            return False
        else:
            return True

    def graphics(self, d2:bool, d3:bool, pd:bool)-> None:

        if d2 == True:
            theta = np.linspace(0, 2*np.pi, 100)
            x = self.center[0] + self.ratio * np.cos(theta)
            y =  self.center[1] + self.ratio * np.sin(theta)

            fig_d2 = plt.figure(figsize=(8,5))
            sns.jointplot(data = self.df[:1][:], x = 0, y = 1)
            plt.axis("on")
            plt.plot(x, y, linestyle = "--", color = "black")
            plt.axhline(self.center[1] - self.radius, linestyle ="--",color = "black")
            plt.axhline(self.center[1] + self.radius, linestyle ="--", color = "black")
            plt.axvline(self.center[0] - self.radius, linestyle ="--", color = "black")
            plt.axvline(self.center[0] + self.radius, linestyle ="--", color = "black")
            plt.text(x = 0, y = 1, s = f"{self.porcentage[0]}")

            plt.title(f"{self.kwargs["d2_title"]}") if self.kwargs["d2_title"] is not None else None

        if d3 == True:

            theta = np.linspace(0, 2*np.pi, 100)
            phi = np.linspace(0, np.pi, 100)

            x = self.center[0] + self.radius * np.cos(theta) * np.sin(phi)
            y = self.center[1] + self.radius * np.sin(theta) * np.sin(phi)
            z = self.center[2] + self.radius * np.cos(phi)

            fig_d3 = plt.figure(figsize = (8, 5))

            ax = fig_d3.add_subplot(projection = "3d")

            ax.scatter(xs = self.df[0][:], ys = self.df[1][:], zs = self.df[2][:], s =1, c = "blue")
            ax.plot3D(xs = x, ys = y, zs = z, linestyle = "--")
            ax.set_title(f"{self.kwargs["d3_title"]}") if self.kwargs["d3_title"] is not None else None

        if pd == True:

            fig_pd = plt.figure(figsize=(8, 5))
            plt.plot(self.porcentage, np.arange(1, self.p + 1, 1), linestyle = "-", linewidth = 2)
            plt.xlabel("Dimension (p)")
            plt.ylabel("Porcentage %")
            plt.title("Porcentage inside the hypersphere vs dimension")

        plt.plot()

class Hypercube():
    """
    Description:
    -------------

    This class calculates the distance of two random samples inside the unit hypercube with a given dimension (p) and number of samples (n). This class works generating a dataset (p x n) and calculating the average distance in each dimension.

    The main purpose of this class is showing both the correlation decreases and the distance increases beetween the same point adding one dimension per iteration. The 

    Parameters:
    ------------

    p: int
        dimension.
    n: int
        number of samples
    correlation: bool
        this parameter controls whether the correlation beetween points is made and plot it or not.
    """

    def __init__(self, p: int, n: int, correlation: bool):

        self.p: int = p
        self.n: int = n
        self.correlation: bool = correlation

    def generate_data(self, p:int, n:int ):

        try: 
            x0 = np.random.uniform(low = 0, high = 1, size = (n, p))
            self.df = pd.DataFrame(x0)
            return self.df
        except:
            raise ValueError("Both the dimension and number of samples must be integers")
      
    def avr_distance(self, p:int, n:int)-> Sequence:

        average_distance = []
        for d in range(p):
            if d == 0: continue
            distance = 0
            for i in n+1:

                if i == n+1: continue
                point_1 = self.df[:d][i]
                point_2 = self.df[:d][i+1]

                distance += self.distance(point1= point_1, point2=point_2)

            average_distance.append(distance)

        return average_distance

    def same_point(self, p:int)->list:

        distance_per_dimension = []

        if self.correlation == True:
            correlation_per_dimension = []

            for d in range(p):
                point_1 = self.df[:d][0]
                point_2 = self.df[:d][1]

                distance_per_dimension.append(self.distance(point1=point_1, point2=point_2))    
                correlation_per_dimension.append(self.correla(point1 = point_1, point2 = point_2))

            return distance_per_dimension, correlation_per_dimension
        
        else:
            for d in range(p):
                point_1 = self.df[:d][0]
                point_2 = self.df[:d][1]
            
                distance_per_dimension.append(self.distance(point1=point_1, point2=point_2))   
            return distance_per_dimension

    def distance(self, point1: Sequence, point2: Sequence) -> float:

        point1 = np.array(point1)
        point2 = np.array(point2)

        distance = np.linalg.norm(point1 - point2)

        return distance

    def correla(self, point1: Sequence, point2:Sequence)->float:    

        point1 = pd.Series([point1])
        point2 = pd.Series([point2])

        cor = point1.corr(point2)

        return cor

    def graphics(self) -> Figure:

        if self.correlation == True:

            fig = plt.figure(figsize = (8,6))
            from matplotlib.gridspec import GridSpec
            gs = GridSpec(2,2)

            avr_d_gr = fig.add_subplot(gs[:, 0])
            corr_gr = fig.add_subplot(gs[0, 1])
            d_gr = fig.add_subplot(gs[1, 1])

            dis, cor = dis =  self.same_point(p = self.p)
            avr_distance = self.avr_distance(p= self.p, n = self.n)

            avr_d_gr.plot(np.arange(1, self.p, 1), avr_distance, linstyle = "-", color = "purple")
            x = np.linspace(0, self.p, 100)
            avr_d_gr.plot(x, np.sqrt(x), linestyle = "--")

            corr_gr.plot(np.arange(1, self.p, 1), cor, linstyle= "-", color = "orange")
            d_gr.plot(np.arange(1, self.p, 1), dis, linestyle = "-", color = "green")

            return fig
        
        else:
            fig_normal = plt.figure(figsize=(8,6))
            ax1 = fig_normal.add_subplot(0, 0)
            ax2 = fig_normal.add_subplot(0, 1)

            avr_dis = self.avr_distance(p =self.p, n = self.n)
            s_point = self.same_point()
            ax1.plot(np.arange(1, self.p + 1, 1), avr_dis, linestyle = "-", color = "purple")
            ax2.plot(np.arange(1, self.p +1, 1), s_point, linestyle = "-", color = "orange")

            return fig_normal

@auditor
class Curse_of_dimensionality():
    """
    Description:
    ------------

    This class provides a reasonable example about how the dimensional scale is counterintuitive.
    The curse of dimensionality, also called, 

    Parameters:
    -----------

    dimension: int
        This parameter controls the dimension space.
    n_sampling: int
        This paramenter controls how many points will generate in the d-dimensional space. A large number will provide better results. See the readme.md
    type_norm: str
        This parameter controls which kind of distance use. Distance available so far: ["Euclidean"]
    graphics: bool
        This parameter controls wheter plot the results in two-dimension graphics: number of sampling vs distance average or not. In case of False, the output will be the final distance average.
    """
    def __init__(self, dimension: int, n_sampling: int,  type_norm: str = "Euclidean", graphics: bool = True):
        pass

    def lenght():
        pass

    def hypercubes():
        pass
    def correlation():
        pass