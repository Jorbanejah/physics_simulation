'''
The following program contains two classes related to the curse of dimensionality.

The first class shows how the proportion of points contained inside a hypersphere
decreases as the dimensionality increases.

The second class shows how, as dimensionality increases:
    1. The average distance between random points increases.
    2. The correlation between random vectors tends to decrease.
    3. The distance between the same two points increases as dimensions are added.

The main purpose of this program is both educational and to illustrate one of
the consequences of the curse of dimensionality in Machine Learning:
dimensionality reduction.

Author: Jorge Orbaneja Huerta
'''

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import time
import sys
import functools
from typing import Sequence
from matplotlib.figure import Figure
from scipy.stats import pearsonr


def auditor(func):
    """
    Decorator that measures and prints the execution time of a function
    and its returned result.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        if not args:
            raise ValueError("The function requires an argument.")

        start = time.time()
        score = func(*args, **kwargs)
        end = time.time()

        duration = (end - start) * 1000

        #print(f"Score: {score}")
        print(f"Duration: {duration:.2f} ms")

        return score

    return wrapper


class Hypersphere():
    """
    This class calculates the percentage of p-dimensional data points
    that fall inside a hypersphere with a given center and radius.

    The data is generated from a normal distribution with a given mean
    and standard deviation.

    Parameters
    ----------
    p : int
        Dimension of the dataset.
    n : int
        Number of samples.
    center : list
        Center of the hypersphere.
    radius : int or float
        Radius of the hypersphere.
    graphic_pd : bool, default=True
        If True, the class generates the 2D/3D graphics together with
        the percentage-versus-dimension graphic.

    Optional parameters
    -------------------
    mean : float
        Mean of the normal distribution. Default is 0.
    std : float
        Standard deviation of the normal distribution. Default is 1.
    title : str
        Title of the main figure.
    """

    def __init__(self, p: int, n: int, center: list, radius: int, graphic_pd: bool = True, **kwargs):

        if not isinstance(p, int) or p < 1:
            raise ValueError("The dimension variable (p) must be a positive integer.")

        if not isinstance(n, int) or n < 1:
            raise ValueError("The number of samples (n) must be a positive integer.")

        if len(center) < p:
            raise ValueError("The center must contain at least p coordinates.")

        if radius <= 0:
            raise ValueError("The radius must be greater than zero.")

        self.p: int = p
        self.n: int = n
        self.center: list = center
        self.radius: int = radius
        self.graphic_pd: bool = graphic_pd

        self.kwargss = kwargs

    def generate_data(self):
        """
        Generate a DataFrame containing n samples of p-dimensional data.

        By default, the data is generated from a standard normal distribution
        with mean 0 and standard deviation 1.
        """

        mean = self.kwargss.get("mean", 0)
        scale = self.kwargss.get("std", 1)

        if not isinstance(mean, (int, float)):
            raise ValueError("The mean must be a number.")

        if not isinstance(scale, (int, float)) or scale <= 0:
            raise ValueError("The standard deviation (std) must be a positive number.")

        x0 = np.random.normal(loc=mean, scale=scale, size=(self.n, self.p))

        self.df = pd.DataFrame(x0)

        return self.df

    @auditor
    def stratified(self, df, p: int) -> list[float]:
        """
        Calculate the percentage of points inside the hypersphere
        for every dimension from 1 to p.
        """

        if p < 1:
            raise ValueError("The dimension must be greater than zero.")

        if df.shape[1] < p:
            raise ValueError("The DataFrame does not contain enough dimensions.")

        self.porcentage = []

        for dimension in range(1, p + 1):

            percent = dimension / p * 100
            sys.stdout.write(f"\rProgress: {percent:.1f}%")
            sys.stdout.flush()

            new_df = df.iloc[:, :dimension]

            percentage = self.porcentage_inside(df=new_df, p=dimension, n=self.n)

            self.porcentage.append(percentage)

        print()

        return self.porcentage

    def porcentage_inside(self, df, p: int, n: int) -> float:
        """
        Calculate the percentage of p-dimensional points that fall
        inside the hypersphere.
        """

        if p < 1:
            raise ValueError("The dimension must be greater than zero.")

        if len(self.center) < p:
            raise ValueError("The center does not contain enough coordinates.")

        center_coordinate = self.center[:p]

        points_inside = 0

        for i in range(n):

            point = df.iloc[i, :p].to_numpy()

            if self.in_or_out(center=center_coordinate, radio=self.radius, point=point):
                points_inside += 1

        percentage = points_inside / n * 100

        return percentage

    def in_or_out(self, center: list, radio: int, point: list) -> bool:
        """
        Determine whether a point is inside the hypersphere.
        """

        center = np.asarray(center)
        point = np.asarray(point)

        distance = np.linalg.norm(center - point)

        return distance <= radio

    def graphics(self, percentage: list) -> Figure:

        if self.graphic_pd:

            fig = plt.figure(figsize=(8, 8))

            from matplotlib.gridspec import GridSpec

            sp = GridSpec(2, 2, fig)

            # Since jointplot cannot be directly embedded into a nested
            # figure, create a scatter plot and two marginal histograms
            # manually.

            sp_joint = sp[0, 1].subgridspec(
                2,
                2,
                width_ratios=[1, 1],
                height_ratios=[1, 1],
                wspace=0.05,
                hspace=0.05, 
            )

            ax_joint = fig.add_subplot(sp_joint[1, 0]) #Bottom-left
            ax_margs_x = fig.add_subplot(sp_joint[0, 0])#Upper-right
            ax_margs_y = fig.add_subplot(sp_joint[1, 1])#Bottom- right
            ax_margs_tl = fig.add_subplot(sp_joint[0, 1])#Upper right
            ax_margs_tl.set_visible(False)

            # The rest of the figures.

            ax_3d = fig.add_subplot(sp[1, 1], projection="3d")
            ax_pd = fig.add_subplot(sp[:, 0])

            # ---------------------------------------------------------
            # 2D GRAPHIC
            # ---------------------------------------------------------

            theta = np.linspace(0, 2 * np.pi, 100)

            x = (self.center[0] + self.radius * np.cos(theta))

            y = (self.center[1] + self.radius * np.sin(theta))

            sns.scatterplot(data= self.df.iloc[:, :2], x = 0, y =1, ax = ax_joint)
            sns.histplot(data = self.df.iloc[:, :2], x =0, ax = ax_margs_x)
            sns.histplot(data = self.df.iloc[:, :2], y =0, ax = ax_margs_y)

            ax_margs_x.set_xticklabels([])
            ax_margs_y.set_yticklabels([])

            ax_joint.plot(x, y, linestyle = "--", color = "black")
            ax_joint.axhline(self.center[1] - self.radius, linestyle ="--",color = "black")
            ax_joint.axhline(self.center[1] + self.radius, linestyle ="--", color = "black")
            ax_joint.axvline(self.center[0] - self.radius, linestyle ="--", color = "black")
            ax_joint.axvline(self.center[0] + self.radius, linestyle ="--", color = "black")
            # Establish the axes.

            max_x = np.max(np.abs(self.df.iloc[:, 0]))
            max_y = np.max(np.abs(self.df.iloc[:, 1]))
            maximus = max(max_x, max_y, self.radius)

            ax_joint.set_xlim([-maximus, maximus])
            ax_joint.set_ylim([-maximus, maximus])

            ax_margs_x.set_title(f"Percentage inside: {self.porcentage[0]:.2f} %")

            ax_joint.set_xlabel(r"$x_1$")
            ax_joint.set_ylabel(r"$x_2$")

            # ---------------------------------------------------------
            # 3D GRAPHIC
            # ---------------------------------------------------------

            resolution = 100
            theta, phi = np.meshgrid( np.linspace(0, np.pi, resolution), np.linspace(0, 2*np.pi, resolution))

            x = self.center[0] + self.radius * np.sin(theta) * np.cos(phi)
            y = self.center[1] + self.radius * np.sin(theta) * np.sin(phi)
            z = self.center[2] + self.radius * np.cos(theta)
           
            ax_3d.scatter(xs = self.df[0][:], ys = self.df[1][:], zs = self.df[2][:], s =1, c = "blue")
            ax_3d.plot3D(xs = x, ys = y, zs = z, linestyle = "--", color = "red")

            # Establish the axes.
            max_x = np.max(np.abs(self.df[0]))
            max_y = np.max(np.abs(self.df[1]))
            max_z = np.max(np.abs(self.df[2]))

            maximus = max(max_x, max_y,max_z,self.radius)

            ax_3d.set_xlim([-maximus, maximus])
            ax_3d.set_ylim([-maximus, maximus])
            ax_3d.set_zlim([-maximus, maximus])

            ax_3d.set_xlabel(r"$x_1$")
            ax_3d.set_ylabel(r"$x_2$")
            ax_3d.set_zlabel(r"$x_3$")

            ax_3d.set_title(f"Percentage inside: {self.porcentage[1]:.2f} %")

            # ---------------------------------------------------------
            # PERCENTAGE VS DIMENSION
            # ---------------------------------------------------------
            percentage_new = [per for per in percentage if per != 0]
            for i in range(2):
                percentage_new.append(0)
            dimensions = np.arange(1, len(percentage_new) +1)
            ax_pd.plot(dimensions, percentage_new, linestyle = "-", linewidth = 2, color = "purple")
            ax_pd.plot(dimensions, percentage_new, "*", linewidth = 2, color = "orange")
            ax_pd.axhline(y = 0, linestyle = "--", color = "black")
            ax_pd.set_xlabel("Dimension (p)")
            ax_pd.set_ylabel("Percentage %")

            title = self.kwargss.get("title")

            if title is not None:
                fig.suptitle(title)

            plt.tight_layout()
            return fig

        else:

            fig_pd = plt.figure(figsize=(8, 5))
            percentage_new = [per for per in percentage if per != 0]
            for i in range(2):
                percentage_new.append(0)

            dimensions = np.arange(1, self.p + 1)
            plt.plot(dimensions, percentage_new, linestyle = "-", linewidth = 2, color = "purple")
            plt.xlabel("Dimension (p)")
            plt.ylabel("Percentage %")

            title = self.kwargss.get("title")

            if title is not None:
                plt.title(title)

            return fig_pd


class Hypercube():
    """
    This class calculates distances between random samples inside
    the unit hypercube as dimensionality increases.

    It shows three effects associated with the curse of dimensionality:

        1. The average distance between random points increases.
        2. The correlation between random vectors tends toward zero.
        3. The distance between the same two points increases as dimensions
           are added.

    Parameters
    ----------
    p : int
        Dimension of the hypercube.
    n : int
        Number of random samples.
    """

    def __init__(self, p: int, n: int):

        if not isinstance(p, int) or p < 1:
            raise ValueError("The dimension variable (p) must be a positive integer.")

        if not isinstance(n, int) or n < 1:
            raise ValueError("The number of samples (n) must be a positive integer.")

        self.p: int = p
        self.n: int = n

    def generate_data(self):
        """
        Generate n random points inside a p-dimensional unit hypercube.
        Every coordinate is sampled independently from Uniform(0, 1).
        """

        x0 = np.random.uniform(low=0,high=1,size=(self.n, self.p))

        self.df = pd.DataFrame(x0)

        return self.df

    @auditor
    def avr_distance(self,df, p: int,n: int) -> Sequence:

        average_distance = []
        correlation = []
        max_dis = []
        min_dis = []

        for d in range(1, p + 1):

            percent = d / p * 100

            sys.stdout.write(f"\rProgress: {percent:.1f}%" )
            sys.stdout.flush()

            distance = []
            corr_values =[]
            for _ in range(n):

                point1, point2 = np.random.choice(a = self.n, size = 2, replace = False)

                point_1 = self.df.iloc[point1,:d]
                point_2 = self.df.iloc[point2,:d]

                idx = np.random.choice(a = self.p, size =2, replace = False)
                vector1 = df.iloc[:d, idx[0]] - df.iloc[:d, idx[0]].mean()
                vector2 = df.iloc[:d, idx[1]] - df.iloc[:d, idx[1]].mean()

                distance.append(self.distance(point1= point_1, point2=point_2))
                cor = self.correla(vector1=vector1, vector2= vector2)

                if cor is not None:
                    corr_values.append(cor)

            max_dis.append(max(distance))
            min_dis.append(min(distance))

            average_distance.append(sum(distance) / n)

            if corr_values:
                correlation.append(sum(corr_values) / len(corr_values))

        print()

        return (average_distance,correlation,max_dis,min_dis)

    @auditor
    def same_point(self, p:int)->list:

        distance_per_dimension = []

        p1, p2  = np.random.sample((2,))

        for d in range(p+1):
            porcentage = d/p*100
            sys.stdout.write(f"\rProgress: {porcentage:.1f} %")
            sys.stdout.flush()

            if d == 0: continue

            point_1 = self.df.iloc[int(self.n * p1), :d]
            point_2 = self.df.iloc[int(self.n * p2), :d]
            
            distance_per_dimension.append(self.distance(point1=point_1, point2=point_2)) 

        print()#new line  
        return distance_per_dimension

    def distance(self, point1: Sequence,point2: Sequence) -> float:

        point1 = np.asarray(point1)
        point2 = np.asarray(point2)

        return np.linalg.norm(point1 - point2)

    def correla( self,vector1: Sequence,vector2: Sequence) -> float:

        if len(vector1) <= 1:
            return None

        cor = (vector1 * vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
        cor = np.abs(cor).mean()
        return cor

    def graphics(self) -> Figure:

        fig = plt.figure(figsize=(8, 6))

        from matplotlib.gridspec import GridSpec

        gs = GridSpec(2, 2)

        avr_d_gr = fig.add_subplot(gs[:, 0])
        corr_gr = fig.add_subplot(gs[0, 1])
        d_gr = fig.add_subplot(gs[1, 1])

        #Functions
        dis =  self.same_point(p = self.p)
        avr_distance, cor, max_dis, min_dis = self.avr_distance(p= self.p, n = self.n)

        dimensions = np.arange(1, self.p + 1)

        # ---------------------------------------------------------
        # AVERAGE DISTANCE
        # ---------------------------------------------------------
        arange_avr_distance = np.arange(1, self.p+1, 1)
        avr_d_gr.set_xscale("log")
        avr_d_gr.plot(arange_avr_distance, avr_distance, linestyle = "-", color = "purple", label = "Avrg")
        avr_d_gr.plot(arange_avr_distance, max_dis, "--", color = "black", label = "max/min")
        avr_d_gr.plot(arange_avr_distance, min_dis, "--", color = "black")
        avr_d_gr.set_xlabel("Dimension (p)")
        avr_d_gr.set_ylabel("Average distance")
        avr_d_gr.set_title("Avrg. distance through dimension")
        avr_d_gr.legend()

  
        # ---------------------------------------------------------
        # CORRELATION
        # ---------------------------------------------------------

        corr_gr.set_xscale("log")

        # Pearson correlation is undefined for one-dimensional
        # vectors, so correlation starts at dimension 2.

        correlation_dimensions = np.arange(2,self.p + 1)

        #Correlation:
        corr_gr.set_xscale("log")
        corr_gr.plot(correlation_dimensions, cor, linestyle= "-", color = "purple")
        corr_gr.set_xlabel("Dimension")
        corr_gr.set_ylabel("E[|r|]")
        corr_gr.set_title("Correlation through dimension")

        # ---------------------------------------------------------
        # DISTANCE BETWEEN THE SAME TWO POINTS
        # ---------------------------------------------------------

        d_gr.plot(arange_avr_distance, dis, linestyle = "-", color = "green")
        d_gr.set_xlabel("Dimension")
        d_gr.set_ylabel("Distance")
        d_gr.set_title("Distance through dimension")

        plt.tight_layout()
        return fig


if __name__ == "__main__":

    dimension = 100
    n_samples = 1000

    Hp = Hypersphere(
        p = dimension,
        n = n_samples,
        center = np.repeat(0, 100),
        radius = 1,
    )

    df = Hp.generate_data()

    percentage = Hp.stratified(df, p = dimension)

    fig = Hp.graphics(percentage=percentage)
    import os
    directory = os.getcwd()
    path = os.path.join(directory, "figures\\hypersphere.png")
    plt.savefig(path,dpi =300, bbox_inches = "tight" )