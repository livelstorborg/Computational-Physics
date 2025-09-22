#include <iostream>
#include <armadillo>
#include <string>
#include <iomanip>
#include <ctime> 
#include <fstream> 
#include "problem2_functions.hpp"
#include "problem3_functions.hpp"
#include "problem4_functions.hpp"
#include "problem5_functions.hpp"

int main() {
    // Seed the random number generator with the current time
    arma::arma_rng::set_seed(std::time(nullptr)); // Set the seed for random number generation

    // a)
    double eps = 1.0e-8;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    const int maxiter = 10000;
    int iterations;
    bool converged;
    arma::vec N = {5, 10, 15};

    int width = 20;

    std::string filename = "N_vs_similarity_transformation.txt";
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::left << std::setw(width) << "N" << std::setw(width) << "similarity transformation" << std::endl;

    for (int i = 0; i < N.size(); ++i) {
        arma::mat A = set_up_A_matrix(N[i]);

        jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

        // Printing
        A.print();
        std::cout << "\n";

        ofile << std::left << std::setw(width) << N[i] << std::setw(width) << iterations << std::endl; // Making header for file
    }

    ofile.close();
    

    return 0;
}
