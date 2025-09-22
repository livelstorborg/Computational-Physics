#ifndef __problem4_functions_hpp__
#define __problem4_functions_hpp__
#include <armadillo>

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int& iterations);

void test_jacobi_rotate(int N);

#endif
