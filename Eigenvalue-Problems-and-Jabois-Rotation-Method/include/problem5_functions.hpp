#ifndef __problem5_functions_hpp__
#define __problem5_functions_hpp__
#include <armadillo>

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);

#endif
