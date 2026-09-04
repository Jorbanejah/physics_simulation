'''
The following program contains a class calls as the file.

This class shows us how high-dimensional space makes bigger when the dimension increase. 

Given the dimension (d), the numbers of point (N) and the norm type, the following class will calculate the distance (norm type) between two random 
d-dimensional points N times and plot it in a two-dimensional graphic (N vs avergare distance between points).

The main purpose of this class are both curiosity and show one performance of Machine Learning: dimensionaly reduction.
'''

import time

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

@auditor
class distance_function(dimension: int, n_sampling: int,  type_norm: str = "Euclidean", graphics: bool = True):
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
    def __init__():
        pass
    pass