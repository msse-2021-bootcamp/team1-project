from datetime import time
import mcsim as mc
import numpy as np
import math
import matplotlib.pyplot as plt
import time


"""Tools for calculating radial distribution functions"""
## You'll need to import your calculate distance function.
def get_all_pairwise_distances(configuration, box_length):
    """
    Get the distances for all the particles.

    Parameters
    ----------
    configuration: np.ndarray
        A set of coordinates
    box_length: float
        The side length of the box

    Returns
    -------
    distances: np.ndarray
        A one dimensional numpy array containing all of particle-particle distances.
    """

    distances = np.array([])
    num_atoms = len(configuration)
    for i in range(num_atoms):
        atom_distances = mc.calculate_distance_np(configuration[i], configuration, box_length)
        atom_distances = np.delete(atom_distances, i)

        distances = np.append(distances, atom_distances)
    return distances


def rdf(values, n_bins, max_value, num_particles, box_length):
    """
    Compute the RDF for a set of values

    Parameters
    -----------
    values : np.ndarray
        The distance values to compute the RDF for.
    n_bins : int
        the number of bins to use for the histogram
    max_value : float
        The maximum value for which to compute the RDF.
    num_particles : int
        The number of particles in the system.
    box_length : float
        The box length

    Returns
    -------
    bin_centers : np.ndarray
        An array of the bin centers
    rdf : np.ndarray
        An array of the rdf values.

    """

    histogram, bins = np.histogram(values, bins=n_bins, range=(0, max_value))

    bin_size = bins[1] - bins[0]

    bin_centers = bins + bin_size / 2
    bin_centers = bin_centers[:-1]

    rdf = []

    rdf = histogram / (
        4 * math.pi * bin_centers ** 2 * bin_size * num_particles ** 2 / box_length ** 3
    )

    return bin_centers, rdf


def time_average_rdf(configurations, box_length, num_atoms, max_value=None):
    """Calculates the average RDF given a set of configurations.
    
    Parameters
    ----------
    configurations : list
        A list containing system coordinates. Can contain multiple frames.
        Arrays should be size (num_particles, 3, num_frames).
    box_length : float
        The box length.
    num_atoms : int
        The number of atoms in the system.
    max_value: float
        The maximum value to consider for the RDF measurement.

    Returns
    -------
    bins: np.n dar
        A list containing the midpoint for each bin.
    avg_rdf: np.ndarray
        The average RDF across all the frames.
    rdfs: np.ndarray
        An array containing the RDF for each frame.
    """

    if not max_value:
        max_value = box_length/2

    rdfs = []
    for configuration in configurations:
        distances = get_all_pairwise_distances(configuration, box_length)

        bins, configuration_rdf = rdf(distances, 200, max_value, num_atoms, box_length)
        rdfs.append(configuration_rdf)

    rdfs = np.array(rdfs)

    avg_rdf = rdfs.mean(axis=0)

    return bins, avg_rdf, rdfs

def main():
    reduced_temp_liquid = .854369559
    reduced_density_liquid = .8772848

    reduced_temp_vapor = .75722972973
    reduced_density_vapor = .0035242

    # reduced_temp_vapor = 3.92209459459
    # reduced_density_vapor = 0.000656010116

    timestamp = time.time()
    print(timestamp)


    cutoff = 3.0

    # starting_coords_liquid, box_length_liquid = mc.read_xyz('../../lj_sample_configurations/lj_sample_config_periodic1.txt')
    starting_coords_liquid, box_length_liquid = mc.generate_config(500,reduced_density_liquid)
    starting_coords_vapor, box_length_vapor = mc.generate_config(500,reduced_density_vapor)

    equilibration_coords_vapor = mc.run_simulation_np(np.array(starting_coords_vapor),box_length_vapor,cutoff,reduced_temp_vapor,500001,.01,20000)
    sampled_coords_vapor = mc.run_simulation_np(np.array(equilibration_coords_vapor[-1]),box_length_vapor,cutoff,reduced_temp_vapor,2000001,.01,10000)

    equilibration_coords_liquid = mc.run_simulation_np(np.array(starting_coords_liquid),box_length_liquid,cutoff,reduced_temp_liquid,500001,.01,50000)
    sampled_coords_liquid = mc.run_simulation_np(np.array(equilibration_coords_liquid[-1]),box_length_liquid,cutoff,reduced_temp_liquid,2000001,.01,10000)

    # equilibration_coords_vapor = mc.run_simulation_np(np.array(starting_coords_vapor),box_length_vapor,cutoff,reduced_temp_vapor,10001,.01,500)
    # sampled_coords_vapor = mc.run_simulation_np(np.array(equilibration_coords_vapor[-1]),box_length_vapor,cutoff,reduced_temp_vapor,10001,.01,100)

    # equilibration_coords_liquid = mc.run_simulation_np(np.array(starting_coords_liquid),box_length_liquid,cutoff,reduced_temp_liquid,10001,.01,500)
    # sampled_coords_liquid = mc.run_simulation_np(np.array(equilibration_coords_liquid[-1]),box_length_liquid,cutoff,reduced_temp_liquid,10001,.01,100)

    bins_liquid, avg_rdf_liquid, rdfs_liquid = time_average_rdf(sampled_coords_liquid,box_length_liquid,500)
    bins_vapor, avg_rdf_vapor, rdfs_vapor = time_average_rdf(sampled_coords_vapor,box_length_vapor,500)

    rdf_fig1 = plt.figure(figsize=(6,6))
    plot1 = rdf_fig1.add_subplot()
    plot1.plot(bins_liquid,avg_rdf_liquid)
    rdf_fig1.savefig(f"./simulation_output/liquid_methane_production_run_{str(timestamp)}.pdf")


    print("saved production liquid simulation results")


    
    rdf_fig2 = plt.figure(figsize=(6,6))
    plot2 = rdf_fig2.add_subplot()
    plot2.plot(bins_vapor,avg_rdf_vapor)
    rdf_fig2.savefig(f"./simulation_output/vapor_methane_production_run_{str(timestamp)}.pdf")
    

    print("saved production vapor simulation results")


    bins_liquid_eq, avg_rdf_liquid_eq, rdfs_liquid_eq = time_average_rdf(equilibration_coords_liquid,box_length_liquid,500)
    bins_vapor_eq, avg_rdf_vapor_eq, rdfs_vapor_eq = time_average_rdf(equilibration_coords_vapor,box_length_vapor,500)

    
    rdf_fig3 = plt.figure(figsize=(6,6))
    plot3 = rdf_fig3.add_subplot()
    plot3.plot(bins_liquid_eq,avg_rdf_liquid_eq)
    rdf_fig3.savefig(f"./simulation_output/liquid_methane_equilibration_run_{str(timestamp)}.pdf")

    rdf_fig4 = plt.figure(figsize=(6,6))
    plot4 = rdf_fig4.add_subplot()
    plot4.plot(bins_vapor_eq,avg_rdf_vapor_eq)
    rdf_fig4.savefig(f"./simulation_output/vapor_methane_equilibration_run_{str(timestamp)}.pdf")

    print("saved equilibration simulation results")

if __name__=="__main__":
    main()