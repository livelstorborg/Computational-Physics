#include <iostream>
#include <armadillo>
#include "problem2_functions.hpp"
#include "problem3_functions.hpp"
#include "problem4_functions.hpp"
#include "problem5_functions.hpp"

// a)
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){

    /* This function is writing the number of similarity transformations required for matrices with 
    different size N to a table*/

    arma::mat R = arma::eye(A.n_rows, A.n_rows); // Initialize eigenvectors to the identity matrix
    arma::mat A_copy = A;

    converged = false;
    int k;
    int l;
    iterations = 0;

    // Iterate up to the maximum number of iterations or until convergence
    for (int i = 0; i < maxiter; i++) {
        double max_offdiag = max_offdiag_symmetric(A_copy, k, l);

        // If the largest off-diagonal element is smaller than eps, set convergance to true
        if (fabs(max_offdiag) < eps) {
            converged = true;
            std::cout << "Size of matrix: " << A.n_rows << "X" << A.n_rows << " rotations: " << iterations << std::endl;
            eigenvalues = A_copy.diag();
            eigenvectors = R;
            break;
        }
        

        jacobi_rotate(A_copy, R, k, l, iterations); //rotate

    }
}