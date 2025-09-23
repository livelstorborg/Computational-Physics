#ifndef __problem3_functions_hpp__
#define __problem3_functions_hpp__
#include <armadillo>

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

void test_max_offdiag_symmetric(int N);

#endif
