#ifndef __IsingModel_hpp__
#define __IsingModel_hpp__

#include <iostream>
#include <armadillo>
#include <vector>

const long double k_b = 1.0; //1.380649e-23;

class IsingModel
{
    public:
    	int L;                  // Lattice dimension
        double T;               // Temperature
        double J;               // Coupling constant
        arma::Mat<int> spins;   // Matrix for spin configuration

        double average_energy;
        double average_magnetisation;
        double average_energy2;
        double average_magnetisation2;
        double specific_heat;
        double susceptibility;


    // Constructor
    IsingModel(int L_in = 10, double temp_in = 1.0, double J_in = 1.0, bool ordered=false);

    // Function to calculate dE
    double delta_energy(int i, int j);

    void monte_carlo_step();

    void metropolis(int num_steps, std::vector<double>& energies, std::vector<double>& cumulative_energies, std::vector<double>& magnetisation);

    // Function that returns the total energy of the system
    double total_energy(const arma::Mat<int>& spin_config);

    // Function that returns the energy per spin
    double energy_per_spin();

    // Function that returns the total magnetisation of the system
    double magnetisation();

    // Function that returns the magnetisation per spin
    double magnetisation_per_spin();

    double partition_function();

    double probability_state();

};

#endif
