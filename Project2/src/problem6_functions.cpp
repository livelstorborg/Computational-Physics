#include <iostream>
#include <armadillo>
#include <vector>
#include <iomanip>
#include <fstream>
#include "problem2_functions.hpp"
#include "problem3_functions.hpp"
#include "problem4_functions.hpp"
#include "problem5_functions.hpp"
#include "problem6_functions.hpp"

void custom_selection_sort(arma::vec& arr, int& first_index, int& second_index, int& third_index) {
    /*
    Function that sorts an array from lowest to highest value. Used for sorting eigenvalue vector.
    */
    int n = arr.size();
    arma::vec original_indices(n);

    for (int i = 0; i < n; ++i) {
        original_indices[i] = i;
    }

    for (int i = 0; i < n - 1; ++i) {
        int minIndex = i;

        for (int j = i + 1; j < n; ++j) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }

        std::swap(arr[i], arr[minIndex]);
        std::swap(original_indices[i], original_indices[minIndex]);

        if (i == 0) {
            first_index = original_indices[i];
        } else if (i == 1) {
            second_index = original_indices[i];
        } else if (i == 2) {
            third_index = original_indices[i];
            break;
        }
    }
}

void write_to_file_numerical(const int N) {
    /*
    Function that writes the eigenvectors with the three smallest eigenvalues to file.
    */
    double eps = 1.0e-8;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    const int maxiter = 10000;
    int iterations;
    bool converged;
    int width = 30;
    int prec = 15;

    arma::mat A = set_up_A_matrix(N);
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    int first_index, second_index, third_index;
    custom_selection_sort(eigenvalues, first_index, second_index, third_index);

    // Write numerical to file
    arma::mat reduced(N, 3); 
    reduced.col(0) = eigenvectors.col(first_index);
    reduced.col(1) = eigenvectors.col(second_index);
    reduced.col(2) = eigenvectors.col(third_index);

    std::ofstream ofile;
    ofile << std::scientific << std::setprecision(prec);
    std::string numerical_filename = "numerical_smallest_lambda.txt";
    ofile.open(numerical_filename);

    ofile << std::left << std::setw(width) << "numerical1" 
          << std::setw(width) << "numerical2" 
          << std::setw(width) << "numerical3" << std::endl;

    for (int i = 0; i < N; ++i) { 
        ofile << std::left << std::setw(width) << std::setprecision(prec) << reduced(i, 0) 
              << std::setw(width) << std::setprecision(prec) << reduced(i, 1) 
              << std::setw(width) << std::setprecision(prec) << reduced(i, 2) 
              << std::endl;
    }
}


