"""
Tests for mcsim package.
"""

import mcsim.monte_carlo as mc
import math
import numpy as np
import random


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

def test_calculate_distance_np_1():

    points1 = np.array([0, 0, 0])
    points2 = np.array([[0, 1, 0],[1,0,0]])

    expected = np.array([1,1])

    observed = mc.calculate_distance_np(points1, points2)

    assert np.isclose(expected.all(),observed.all())

def test_calculate_distance_np_period_boundaries():

    points1 = np.array([0, 8, 0])
    points2 = np.array([[0, 0, 0]])

    box_length = 10

    expected = 2
    
    observed = mc.calculate_distance_np(points1, points2, box_length)

    assert math.isclose(expected,observed)

def test_accept_or_reject_false():
    delta_energy = 1
    beta = 1
    random.seed(0)

    accepted = mc.accept_or_reject(delta_energy, beta)

    random.seed()

    assert accepted is False

def test_accept_or_reject_true():
    delta_energy = 1
    beta = 1
    random.seed(1)

    accepted = mc.accept_or_reject(delta_energy, beta)

    random.seed()

    assert accepted is True

def test_calculate_pair_energy():
    test_coords = [[0,0,0], [0, 0, 2**(1/6)], [0 ,0 , 2*2**(1/6)]]
    assert mc.calculate_pair_energy(test_coords, 1, 10, 3) == -2

def test_calculate_total_energy():
    coordinates = [[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]]
    assert math.isclose(mc.calculate_total_energy(coordinates,None,3.0), -2, rel_tol=0.05)

def test_calculate_total_energy_np():
    coordinates = np.array([[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]])
    assert np.isclose(mc.calculate_total_energy(coordinates,None,3.0), -2, rtol=0.05)

def test_calculate_LJ():
    assert mc.calculate_LJ(math.pow(2, (1/6))) == -1

def test_calculate_LJ_np():
    assert mc.calculate_LJ_np(np.array([math.pow(2, (1/6)),-1])).all() == np.array([-1,0]).all()

def test_calculate_tail_correction():
    num_particles = 10.0
    box_length = 1.0
    cutoff = 1.0

    expected = (8.0*math.pi/3.0)*100.0*(-2.0/3.0)
    assert math.isclose(mc.calculate_tail_correction(num_particles,box_length,cutoff),expected)

def test_calculate_pair_energy():
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]
    assert mc.calculate_pair_energy(coordinates, 1, 10.0, 3.0) == -2


def test_calculate_pair_energy_np():
    coordinates = np.array([[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]])
    assert mc.calculate_pair_energy_np(coordinates, 1, 10.0, 3.0) == -2

