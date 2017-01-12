from helper import smallest_distance, rotation
import numpy as np
import math


def test_smallest_distance_between_two_points():
	assert smallest_distance([(1, 0, 0)], [(0, 0, 0)]) == 1

def test_smallest_distance_between_small_sets():
	assert smallest_distance([(1, 0, 0)], [(0, 0, 0), (-1, 0, 0)]) == 1
	assert smallest_distance([(0, 0, 0)], [(0, 0, 0), (-1, 0, 0)]) == 0

	assert smallest_distance([(0, 0, 0), (-1, 0, 0)], [(1, 0, 0)]) == 1
	assert smallest_distance([(0, 0, 0), (-1, 0, 0)], [(0, 0, 0)]) == 0


def test_rotate():
	identity = ([1,0,0],[0,1,0],[0,0,1])
	assert np.all(rotation([0,0,0]) == identity)

def test_rotate_z90():
	z90 = ([-1,0,0],[0,-1,0],[0,0,1])
	print(z90)
	answer = np.round((rotation([0,0,math.pi])),decimals=0)
	print(answer)
	assert np.all(answer == z90)