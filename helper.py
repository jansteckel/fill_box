import numpy as np
from math import sqrt, cos, sin

# define the function for the rotation matrix
def rotation(theta):
    tx,ty,tz = theta
    Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
    Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
    Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])
    return np.dot(Rx, np.dot(Ry, Rz))

def smallest_distance(array_a, array_b):
    """
    define a function that evaluates 
    the smallest distance between any two points in two arrays 
    only works for arrays of shape n by 3
    """
    smallest_d = 1000 # big arbitrary start value    
    for a in array_a: 
        for b in array_b:
            dist = sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2) 
            smallest_d = min(smallest_d, dist)
            return smallest_d