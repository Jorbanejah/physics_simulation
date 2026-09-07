'''
The following program contains a class calls as the file.

This class shows us how high-dimensional space makes bigger when the dimension increase. 

Given the dimension (d), the numbers of point (N) and the norm type, the following class will calculate the distance (norm type) between two random 
d-dimensional points N times and plot it in a two-dimensional graphic (N vs avergare distance between points).

The main purpose of this class are both curiosity and show one performance of Machine Learning: dimensionaly reduction.
'''

import time
import numpy as np
from typing import Sequence

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

def ecludian_distant(dimension:int, point_1:Sequence, point_2:Sequence)->float:
    "Calculus of the ecludian distance of two given points in a p-dimension space"

    Xp = 0
    for p in dimension-1:
        Xp += (point_1[p] - point_2[p])**2

    return np.sqrt(Xp)

def hypercircle(dimension: int, radium: float, center: Sequence, point: Sequence)-> bool:

    Xp = ecludian_distant(dimension=dimension, point_1 = center, point_2 = point)

    if Xp >= radium:
        return False
    else:
        return True

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


class Hypercircles():
    """
    Description:
    -------------
    The class calculates the percent of samples that they are within the hypercircle of center(X) and radius r.

    The class performance with two outside function: ecludian_distant that calculates the distance between the sample and the center, and the hypercircle that create it 
    and its ouput tells whether the points is inside or not. The class works splitting the dataset in diferent subset of lower dimensionality, then describing the hypercircle and  
    look for the points that are within it.

    The main purpose of this class is to verify that data tend to be sparse in higher dimensions. So the output, 

    Parameters:
    -------------
    dataset:Sequence
        The dataset contains the dimension and any sample of it.
    Center: Sequence
        A vector of float numbers that describes the center of the hypercircle
    Radius: float
    """
    def __init__(self, dataset: Sequence, center: Sequence, radius: float):
        pass
    def graphics():
        pass
