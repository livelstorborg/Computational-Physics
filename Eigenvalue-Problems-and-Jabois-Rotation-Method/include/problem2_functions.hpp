#ifndef __problem2_functions_hpp__
#define __problem2_functions_hpp__
#include <armadillo>

arma::mat set_up_A_matrix(const int N);

arma::vec solve_eigval(arma::mat A);

arma::mat solve_eigvec(arma::mat A);

void analytical_eig_vec_val(arma::mat& analytical_eigvec, arma::vec& analytical_eigval, int N);

void test_eigval_eigvec(int N);

#endif
