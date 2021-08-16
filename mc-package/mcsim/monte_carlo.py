import math
import random
import numpy as np

"""
Functions For running a Monte Carlo Simulation
"""

def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Calculate the total energy of a set of particles using the Lennard Jones potential.
    
    Parameters
    ----------
    coordinates : list
        A nested list containing the x, y,z coordinate for each particle
    box_length : float
        The length of the box. Assumes cubic box.
    cutoff : float
        The cutoff length
    
    Returns
    -------
    total_energy : float
        The total energy of the set of coordinates.
    """
    
    total_energy = 0
    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            # Calculate the distance between the particles - exercise.
            dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length)

            if dist_ij < cutoff:
                # Calculate the pairwise LJ energy
                LJ_ij = calculate_LJ(dist_ij)

                # Add to total energy.
                total_energy += LJ_ij
    return total_energy


def calculate_total_energy_np(coordinates, box_length, cutoff):
    """
    Calculate the total energy of a set of particles using the Lennard Jones potential.
    
    Parameters
    ----------
    coordinates : list
        A nested list containing the x, y,z coordinate for each particle
    box_length : float
        The length of the box. Assumes cubic box.
    cutoff : float
        The cutoff length
    
    Returns
    -------
    total_energy : float
        The total energy of the set of coordinates.
    """
    
    num_atoms = len(coordinates)

    pair_energies = np.array([calculate_pair_energy_np(coordinates,i,box_length,cutoff) for i in range(num_atoms)])
    
    return pair_energies.sum()


def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.
    
    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.
       
    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates
    """
    
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


def read_xyz_np(filepath):
    """
    Reads coordinates from an xyz file using functions from the python standard library.
    
    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.
       
    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates
    """
    coordinates = np.genfromtxt(filepath,skip_header=2,usecols=[1,2,3])

    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())

    return coordinates, box_length


def calculate_LJ(r_ij):
    """
    The LJ interaction energy between two particles.

    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.
    
    Returns
    -------
    pairwise_energy : float
        The pairwise Lennard Jones interaction energy in reduced units.

    Examples
    --------
    >>> calculate_LJ(1)
    0

    """
    
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    
    return pairwise_energy


def calculate_LJ_np(r_ij):
    """
    The LJ interaction energy between two particles.
    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.
    
    Returns
    -------
    pairwise_energy : float
        The pairwise Lennard Jones interaction energy in reduced units.

    Examples
    --------
    >>> calculate_LJ(1)
    0

    """
    
    r6_term = np.power(1/r_ij, 6)
    r12_term = np.power(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    return pairwise_energy

    
def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two points. When `box_length` is set, the minimum image convention is used to calculate the distance between the points.

    Parameters
    ----------
    coord1, coord2 : list
        The coordinates of the points, [x, y, z]
    
    box_length : float, optional
        The box length. This function assumes box is a cube.

    Returns
    -------
    distance : float
        The distance between the two points.
    """
    
    distance = 0
    
    for i in range(3):
        dim_dist = (coord1[i] - coord2[i])
        if box_length:
            dim_dist = dim_dist - box_length * round(dim_dist/ box_length)
        dim_dist = dim_dist**2
        distance += dim_dist
        
        
    distance = math.sqrt(distance)
    
    return distance


def calculate_distance_np(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.
    Parameters
    ----------
    coord1, coord2: np.array
        The atomic coordinates
    
    box_length : float
        The box length. If given, the minimum image convention will be used to calculate the distance.
    Returns
    -------
    distance: float
        The distance between the two points.
    """
    
    coord_dist = coord1 - coord2
    
    if box_length:
            coord_dist = coord_dist - box_length * np.round(coord_dist / box_length)

    #how many axis it has 
    if coord_dist.ndim < 2:
        coord_dist = coord_dist.reshape(1,-1)

    coord_dist = coord_dist ** 2
    coord_dist_sum = coord_dist.sum(axis=1)
    distance = np.sqrt(coord_dist_sum)
    return distance


def calculate_tail_correction(n,b,r_c):
    r3_term = math.pow(1/r_c, 3)
    r9_term = (1/3) * (math.pow(r3_term, 9))
    NV_term = (8*math.pi/3) * ((n**2)/(b**3))
    tail_correction = NV_term * (r9_term - r3_term)
    
    return tail_correction


def accept_or_reject(delta_U, beta):
    """
    Accept or reject a move based on metropolis criterion.
    
    Paramaters
    ----------------
    delta_U : float
    
        change in energy for moving systems from m to n

    beta : float
        1/temperature

    Returns
    --------
    bool
        whether move is accepted-true

    """
    
    
    if delta_U <= 0.0:
        accept = True
    else:
        #Gen random number on(0,1)
        random_number = random.random()
        p_acc = math.exp(-beta*delta_U)
        
        if random_number < p_acc:
            accept = True
         
        else:
            accept = False
            
    return accept


def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    
    """
    Calculate interaction energy of particle w/ its environment (all other particles in sys)
    
    Parameters
    ----------------
    coordinates : list
        the coordinates for all particles in sys
    i_particle : int
        particle number for which to calculate energy
        
    cutoff : float
        simulation cutoff. beyond distances, interactions aren't calculated
    
    box length : float
    
        length of simultion box. assumes cubic boc
        
    Returns
    ---------------
    
    float
        pairwise interaction energy of ith particle w/all other particles in sys

    """
    
    e_total = 0.0
    i_position = coordinates[i_particle]
    
    num_atoms = len(coordinates)
    
    for j_particle in range(num_atoms):
        if i_particle != j_particle:
            
            j_position = coordinates[j_particle]
            rij = calculate_distance(i_position, j_position, box_length)
            
            if rij < cutoff:
                e_pair = calculate_LJ(rij)
                e_total += e_pair
                
                

    return e_total


def calculate_pair_energy_np(coordinates, i_particle, box_length, cutoff):
    """
    Calculate interaction energy of particle w/ its environment (all other particles in sys)
    
    Parameters
    ----------------
    coordinates : list
        the coordinates for all particles in sys
    i_particle : int
        particle number for which to calculate energy    
    cutoff : float
        simulation cutoff. beyond distances, interactions aren't calculated 
    box length : float
        length of simultion box. assumes cubic boc
        
    Returns
    ---------------
    float
        pairwise interaction energy of ith particle w/all other particles in sys
        
    """

    particle_distances  = calculate_distance_np(coordinates[i_particle], coordinates[i_particle+1:], box_length) 
    particle_distances_filtered = particle_distances[particle_distances < cutoff]
    return calculate_LJ_np(particle_distances_filtered).sum()


def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000):
    """
    Run a Monte Carlo Simulation with the specified parameters
    
    Parameters
    ----------------
    coordinates : list
        the coordinates for all particles in sys
    box length : float
    
        length of simultion box. assumes cubic boc
    cutoff : float
        simulation cutoff. beyond distances, interactions aren't calculated
    reduced_temperature : float
    num_steps : int
        the number of steps to run the simulation for
    max_displacement: float
    freq : int

    """
    # calculated quantities
    beta = 1/reduced_temperature
    num_particles = len(coordinates)

    # Energy calculations
    total_energy = calculate_total_energy(coordinates, box_length,cutoff)

    tail_correction = calculate_tail_correction(num_particles, box_length, cutoff)

    total_energy += tail_correction

    for step in range(num_steps):
            
        # 1. Randomly pick one of particles
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate in interaction energy of selected particle w/system and store this value
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
    
        # 3. Generate random x,y,z displacement
        
        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)


        #4. Modify coordinate of Nth particle by genrated displacements
        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand
        
        #5. Calculate interaction energy of moved particle w/system
        proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        delta_energy = proposed_energy - current_energy
        #6. Calculate if we accept move based on energy diff
        
        accept = accept_or_reject(delta_energy, beta)
        #7. if accept, move particle
        if accept: 
            total_energy += delta_energy
        else:
            # Move is not accepted, roll back coordinates
            coordinates[random_particle][0] -= x_rand
            coordinates[random_particle][1] -= y_rand
            coordinates[random_particle][2] -= z_rand


        #8. print energy if step is multiple of freq
        
        if step % freq == 0:
            print(step, total_energy/num_particles)

    return coordinates


def run_simulation_np(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement, freq=1000):
    """
    Run a Monte Carlo simulation with the specified parameters. 
    """

    # Calculated quantities

    beta = 1 / reduced_temperature

    num_particles = len(coordinates)

    # Energy calculations 

    total_energy = calculate_total_energy_np(coordinates, box_length, cutoff)

    tail_correction = calculate_tail_correction(num_particles, box_length, cutoff)

    total_energy += tail_correction

    for step in range(num_steps):
            
        # 1. Randomly pick one of particles
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate in interaction energy of selected particle w/system and store this value
        current_energy = calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)
    
        # 3. Generate random x,y,z displacement
        
        xyz_rand = np.random.uniform(-max_displacement, max_displacement, [1,3])

        #4. Modify coordinate of Nth particle by genrated displacements
        coordinates[random_particle] = coordinates[random_particle] + xyz_rand

        #5. Calculate interaction energy of moved particle w/system
        proposed_energy = calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)
        
        delta_energy = proposed_energy - current_energy
        
        #6. Calculate if we accept move based on energy diff       
        accept = accept_or_reject(delta_energy, beta)
        #7. if accept, move particle
        if accept: 
            total_energy += delta_energy
        else:
            # Move is not accepted, roll back coordinates
            coordinates[random_particle] = coordinates[random_particle] - xyz_rand

        #8. print energy if step is multiple of freq
        if step % freq == 0:
            print(step, total_energy/num_particles)
