#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <omp.h>
#include "IsingModel.hpp"
#include <iomanip>  // For std::put_time
#include <ctime>    // For std::localtime
#include <chrono> // For timing



void print_current_time(const std::string& message) 
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::cout << message << ": " << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S") << std::endl;
}





void write_to_file_8(std::string filename, arma::vec temperature, std::vector<double> eps, std::vector<double> eps2, std::vector<double> mag, std::vector<double> mag2, std::vector<double> heat_cap, std::vector<double> sus)
{
    int width = 15;  
    int prec = 8;     

    std::ofstream outfile(filename);

    outfile << std::left << std::setw(width) << "Temperature"
            << std::setw(width) << "Energy"
            << std::setw(width) << "Energy2"
            << std::setw(width) << "Magnetization"
            << std::setw(width) << "Magnetization2"
            << std::setw(width) << "Heat Capacity"
            << std::setw(width) << "Susceptibility" << "\n";

    for (size_t i = 0; i < temperature.n_elem; ++i)
    {
        outfile << std::left << std::setw(width) << std::setprecision(prec) << temperature[i]
                << std::setw(width) << std::setprecision(prec) << eps[i]
                << std::setw(width) << std::setprecision(prec) << eps2[i]
                << std::setw(width) << std::setprecision(prec) << mag[i] 
                << std::setw(width) << std::setprecision(prec) << mag2[i]
                << std::setw(width) << std::setprecision(prec) << heat_cap[i]
                << std::setw(width) << std::setprecision(prec) << sus[i] << "\n";
    }
    outfile.close();
}





int main()
{

    int mc_cycles = 100000;
    double dt = 0.01;
    double J = 1.0;


    
    // --------------- Problem 8 (L = {40, 60, 80, 100}) ---------------
    std::vector<int> lattice_sizes = {40, 60, 80, 100};
    arma::vec temperatures8 = arma::regspace(2.1, dt, 2.4 + dt);

    std::cout << "MONTE CARLO CYCLES: " << mc_cycles << std::endl;
    std::cout << "\n";

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int i = 0; i < lattice_sizes.size(); i++)
    {


        std::vector<double> av_energy;
        std::vector<double> av_energy2;
        std::vector<double> av_magnetisation;
        std::vector<double> av_magnetisation2;
        std::vector<double> sp_heat;
        std::vector<double> sus;

        for (int j = 0; j < temperatures8.n_elem; j++)
        {
            std::cout << std::left << "L=" << lattice_sizes[i] << std::setw(10) << "Temperature: " << temperatures8[j] << std::endl;

            IsingModel model_many(lattice_sizes[i], temperatures8[j], J, false);
            std::vector<double> energies;
            std::vector<double> cumulative_energies;
            std::vector<double> magnetisations;

            model_many.metropolis(mc_cycles, energies, cumulative_energies, magnetisations);

            av_energy.push_back(model_many.average_energy);
            av_energy2.push_back(model_many.average_energy2);
            av_magnetisation.push_back(model_many.average_magnetisation);
            av_magnetisation2.push_back(model_many.average_magnetisation2);
            sp_heat.push_back(model_many.specific_heat);
            sus.push_back(model_many.susceptibility);
        }

        std::string filename = "L" + std::to_string(lattice_sizes[i]) + "_func_of_temp.txt";
        write_to_file_8(filename, temperatures8, av_energy, av_energy2, av_magnetisation, av_magnetisation2, sp_heat, sus);

        print_current_time("Finished lattice size L=" + std::to_string(lattice_sizes[i]));
    }
    

    

    /*
    // --------------- Problem 7 - Timing ---------------
    // ----- Parallel -----
    std::vector<int> lattice_sizes = {4, 6, 8, 10};
    arma::vec temperatures7 = arma::regspace(2.1, dt, 2.4 + dt);

    std::cout << "MONTE CARLO CYCLES: " << mc_cycles << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "\n";

    auto start = std::chrono::high_resolution_clock::now();

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int i = 0; i < lattice_sizes.size(); i++)
    {


        std::vector<double> av_energy;
        std::vector<double> av_energy2;
        std::vector<double> av_magnetisation;
        std::vector<double> av_magnetisation2;
        std::vector<double> sp_heat;
        std::vector<double> sus;

        for (int j = 0; j < temperatures7.n_elem; j++)
        {
            std::cout << std::left << "L=" << lattice_sizes[i] << std::setw(10) << "Temperature: " << temperatures7[j] << std::endl;

            IsingModel model_many(lattice_sizes[i], temperatures7[j], J, false);
            std::vector<double> energies;
            std::vector<double> cumulative_energies;
            std::vector<double> magnetisations;

            model_many.metropolis(mc_cycles, energies, cumulative_energies, magnetisations);

            av_energy.push_back(model_many.average_energy);
            av_energy2.push_back(model_many.average_energy2);
            av_magnetisation.push_back(model_many.average_magnetisation);
            av_magnetisation2.push_back(model_many.average_magnetisation2);
            sp_heat.push_back(model_many.specific_heat);
            sus.push_back(model_many.susceptibility);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time: " << duration.count() << " ms" << std::endl;
    */

    return 0;
}