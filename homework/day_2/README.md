# Day 2 Homework

import numpy as np
import math
import os

def calculate_LJ(r_ij):
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    
    return pairwise_energy
def read_xyz(filepath):
    
    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()
    
    atomic_coordinates = []
    
    for atom in coordinates:
        split_atoms = atom.split()
        
        float_coords = []
        
        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))
            
        atomic_coordinates.append(float_coords)
        
    
    return atomic_coordinates, box_length

config1_file = 'C:\\Users\\91095\\Desktop\\xyz file\\lj_sample_config_periodic1(1).txt'

sample_coords, box_length = read_xyz(config1_file)
def calculate_total_pair_energy(coordinates, cutoff):
    total_energy = 0
    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            # Calculate the distance between the particles - exercise.
            dist_ij = calculate_distance(coordinates[i], coordinates[j])

            if dist_ij < cutoff:
                # Calculate the pairwise LJ energy
                LJ_ij = calculate_LJ(dist_ij)

                # Add to total energy.
                total_energy += LJ_ij
    return total_energy

#change cutoff
print(calculate_total_pair_energy(sample_coords, 0.9))
print(calculate_total_pair_energy(sample_coords, 1))
print(calculate_total_pair_energy(sample_coords, 1.1))
print(calculate_total_pair_energy(sample_coords, 1.2))
