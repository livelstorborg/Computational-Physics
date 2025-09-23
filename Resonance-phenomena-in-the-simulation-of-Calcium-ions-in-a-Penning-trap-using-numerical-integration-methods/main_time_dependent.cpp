#include <iostream>
#include <armadillo>
#include <iomanip>
#include <string>
#include <fstream>
#include "Particle.hpp"
#include "PenningTrap.hpp"

void write_to_file(std::string filename, arma::vec x, arma::vec y)
{
    int width = 30;
    int prec = 15;

    std::ofstream ofile;
    ofile.open(filename);

    ofile << std::scientific << std::setprecision(prec);
    ofile << std::left << std::setw(width) << "x"
                       << std::setw(width) << "y" << std::endl;

    for(int i = 0; i < x.n_elem; i++) 
    {
        ofile << std::left << std::setw(width) << x[i]
                           << std::setw(width) << y[i] << std::endl;
    }

    ofile.close();
}

int main() 
{	
	double mass = 40.078;
	int charge = 1;
	int t_max = 500;
	int n_step = 40000;
	double dt = 500. / n_step;
/*
	arma::vec f = {0.1, 0.4, 0.7};  
	arma::vec omega_v = arma::regspace(0.2, 0.02, 2.5);    
    //omega_v = arma::regspace(0.2, 0.1, 2.5); //Example to run fast

	arma::vec particles_inside_trap(omega_v.n_elem);

	for(int i = 0; i < f.n_elem; i++) 
    {
	    for(int j = 0; j < omega_v.n_elem; j++) 
        {
            std::cout << "omega_iteration: " << j << std::endl;
	        PenningTrap trap = PenningTrap(T, 0.025 * V, 500, false, true, f[i], omega_v[j]);
	        trap.add_random_particle(100, charge, mass);

	        for (int k = 0; k < n_step; k++) 
            {
	            trap.change_time(k * dt);
	            trap.evolve_RK4(dt);
            }
	        
	        particles_inside_trap[j] = trap.count_particles();
	    }

	    std::string filename_particles_inside = "txt_files/f" + std::to_string(f[i]) + ".txt";
	    write_to_file(filename_particles_inside, particles_inside_trap, omega_v);  // Ensure format is correct
	}
    
    // ---------- fine grained - without columb interactions ----------
*/
    arma::vec f_fine = {0.1};  
    arma::vec omega_v_fine = arma::regspace(1.1, 0.002, 1.6);   
    //omega_v_fine = arma::regspace(1.1, 0.01, 1.7);  //Example to run fast

    arma::vec particles_inside_trap_fine(omega_v_fine.n_elem);

    for(int i = 0; i < f_fine.n_elem; i++)
    {
        for(int j = 0; j < omega_v_fine.n_elem; j++) 
        {
            std::cout << "omega_iteration: " << j << std::endl;
            PenningTrap trap_fine = PenningTrap(T, 0.025 * V, 500, false, true, f_fine[i], omega_v_fine[j]);
            trap_fine.add_random_particle(100, charge, mass);

            for(int k = 0; k < n_step; k++) 
            {
                trap_fine.change_time(k * dt);
                trap_fine.evolve_RK4(dt);
            }

            particles_inside_trap_fine[j] = trap_fine.count_particles();
        }

        std::string filename_particles_inside = "txt_files/f_fine_non_interaction" + std::to_string(f_fine[i]) + ".txt";
        write_to_file(filename_particles_inside, particles_inside_trap_fine, omega_v_fine);
    }
    
    // ---------- fine grained - with columb interactions ----------

/*
    arma::vec particles_inside_trap_fine_interactions(omega_v_fine.n_elem);

    for(int i = 0; i < f_fine.n_elem; i++)
    {
        for(int j = 0; j < omega_v_fine.n_elem; j++)
        {
            std::cout << "omega_iteration: " << j << std::endl;
            std::cout << "interaction " << j << std::endl;

            PenningTrap trap_fine_interactions = PenningTrap(T, 0.025 * V, 500, true, true, f[i], omega_v[j]);
            trap_fine_interactions.add_random_particle(100, charge, mass);

            for(int k = 0; k < n_step; k++) 
            {
                trap_fine_interactions.change_time(k * dt);
                trap_fine_interactions.evolve_RK4(dt);
            }

            particles_inside_trap_fine_interactions[j] = trap_fine_interactions.count_particles();
        }

        std::string filename_particles_inside = "txt_files/f_fine_interaction" + std::to_string(f_fine[i]) + ".txt";
        write_to_file(filename_particles_inside, particles_inside_trap_fine_interactions, omega_v_fine);
    }
*/
    return 0;
}

