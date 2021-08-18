#define _USE_MATH_DEFINES

#include <random> // for random numbers
#include <chrono> // for generating the seed
#include <fstream> // for reading/writing files
#include <array>   // for std::array
#include <vector>  // for std::vector
#include <utility> // for std::pair
#include <cmath> 
#include <iostream>

// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

/*! Generate a random double within a given range */
double random_double(double lower_bound, double upper_bound)
{
   std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
   return dist(re);
}

/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double random_integer(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
} 

// Make some types more convenient
typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

std::pair<Coordinates, double> read_xyz(std::string file_path)
{
    // Opens up a file stream for input
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if(!infile.is_open())
    {
        throw std::runtime_error("File path in read_xyz does not exist!");
    }

    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;

    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;

    // now the number of atoms
    infile >> num_atoms;

    // Uncomment to help troubleshoot
    //std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;

    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;

    for(int i = 0; i < num_atoms; i++)
    {
        AtomCoord coord;

        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];

        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}

double calculate_LJ(double r_ij)
{
	double r6_term = std::pow((1/r_ij),6);
	double r12_term = std::pow(r6_term,2);
	double pairwise_energy = 4 * (r12_term - r6_term);

	return pairwise_energy;
}

double calculate_tail_correction(int num_particles, double box_length, double cutoff)
{
	double r3_term = std::pow((1/cutoff),3);
	double r9_term = (1/3) * std::pow(r3_term, 9);
	double NV_term = (8*M_PI/3) * (double)(std::pow(num_particles,2)/std::pow(box_length,3));
	double tail_correction = NV_term * (r9_term - r3_term);

	return tail_correction;
}

bool accept_or_reject(double delta_u, double beta)
{
	bool accept;
	if(delta_u <= 0.0)
	{
		accept = true;
	}
	else
	{
		double random_number = random_double(0.0,1.0);
		double p_acc = std::pow(M_E,(-beta * delta_u));

		if(random_number < p_acc)
		{
			accept = true;
		}
		else
		{
			accept = false;
		}
	}

	return accept;

}

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length)
{
	double distance = 0.0;

	for(int i = 0; i < 3; i++){
		double dim_dist = coord1[i] - coord2[i];
		if(box_length != 0.0){
			dim_dist = dim_dist - box_length * std::round(dim_dist/box_length);
		}
		dim_dist = std::pow(dim_dist,2);	
		distance += dim_dist;
	}
		
	distance = std::sqrt(distance);

	return distance;
}

double calculate_pair_energy(Coordinates coords,int i_particle,double box_length,double cutoff)
{
	double e_total = 0.0;
	AtomCoord i_position = coords[i_particle];
	
	int num_atoms = coords.size();

	for(int j_particle=0; j_particle < num_atoms; j_particle++)
	{ 
		if(i_particle != j_particle)
		{
			AtomCoord j_position = coords[j_particle];
		       	double rij = calculate_distance(i_position,j_position,box_length);

			if(rij < cutoff)
			{
				double e_pair = calculate_LJ(rij);
				e_total += e_pair;
			}
		}
	}

	return e_total;

}

double calculate_total_energy(Coordinates coords, double box_length, double cutoff)
{
	double total_energy = 0.0;

	int num_atoms = coords.size();

	for(int i = 0; i < num_atoms; i++)
	{
		AtomCoord i_position = coords[i];
	        for(int j=i+1; j < num_atoms; j++)
        	{
                       	AtomCoord j_position = coords[j];
                       	double dist_ij = calculate_distance(i_position,j_position,box_length);

                       	if(dist_ij < cutoff)
                 	{
                               	total_energy += calculate_LJ(dist_ij);
                       	}		
       	 	}
		
	}
	return total_energy;
}	

Coordinates run_simulation(Coordinates coords, double box_length, double cutoff, double reduced_temperature, int num_steps, double max_displacement, int freq)
{	
	// calculated quantities
	double beta = 1/reduced_temperature;
	int num_particles = coords.size();

	// energy calculations 
	double total_energy = calculate_total_energy(coords,box_length,cutoff);
	double tail_correction = calculate_tail_correction(num_particles, box_length, cutoff);
	total_energy += tail_correction;

	for(int j = 0; j < num_steps; j++)
	{
		// randomly pick a particle
		int random_particle = random_integer(0,num_particles);

		// calculate interaction energy of selected particles with the system 
		double current_energy = calculate_pair_energy(coords,random_particle,box_length, cutoff);

		//generate random x,y,z displacement
		double x_rand = random_double(-max_displacement,max_displacement);
		double y_rand = random_double(-max_displacement,max_displacement);
		double z_rand = random_double(-max_displacement,max_displacement);

		// modify coordinates of the nth particle by generated displacements
		AtomCoord new_coord;
		AtomCoord old_coord = coords[random_particle];
		new_coord[0] = coords[random_particle][0] + x_rand;
		new_coord[1] = coords[random_particle][1] + y_rand;
		new_coord[2] = coords[random_particle][2] + z_rand;
		coords[random_particle] = new_coord;

		//calculate interaction energy of moved particle with system
		double proposed_energy = calculate_pair_energy(coords, random_particle, box_length, cutoff);
		double delta_energy = proposed_energy - current_energy;

		// calculate if we accept of reject move based on energy diff
		bool accept = accept_or_reject(delta_energy,beta);
		if(accept == true)
		{
			total_energy += delta_energy;
		}
		else
		{
			//move is not accepted roll back change
			coords[random_particle] = old_coord;
		}

		// print energy  if step is multiple of freq

		if(j%freq == 0)
		{
			std::cout << "Energy at Step " << j << ": " << total_energy/num_particles << std::endl;
		}		

	}
	return coords;

}	

int main(void)
{
       	// Initialize random number generation based on time
	re.seed(std::chrono::system_clock::now().time_since_epoch().count());

   	 /* other code here */
   	 std::pair<Coordinates, double> xyz_info = read_xyz("../../lj_sample_configurations/lj_sample_config_periodic1.txt");

   	 Coordinates coords = xyz_info.first;
   	 double box_length = xyz_info.second;

  	 // can now use box_length and coords
	 
	 double test_total = calculate_total_energy(coords,box_length,3.0);

	 std::cout << test_total << std::endl;
	 //set simulation params
	 double reduced_temperature = 1.5;
	 int num_steps = 5000;
	 double max_displacement = 0.1;
	 double cutoff = 3.0;
	 int freq = 1000;

	 //run simulation
	 Coordinates out_coords = run_simulation(coords,box_length,cutoff,reduced_temperature,num_steps,max_displacement,freq);
}




