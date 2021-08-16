"""
Tests for mcsim package.
"""

import mcsim.monte_carlo as mc
import math

def test_calculate_distance_1():

    points1 = [0, 0, 0]
    points2 = [0, 1, 0]

    expected = 1

    observed = mc.calculate_distance(points1, points2)

    assert math.isclose(expected,observed)

def test_calculate_distance_period_boundaries():

    points1 = [0, 8, 0]
    points2 = [0, 0, 0]

    box_length = 10

    expected = 2
    
    observed = mc.calculate_distance(points1, points2, box_length)

    assert math.isclose(expected,observed)