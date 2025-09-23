#ifndef __problem7_functions_hpp__
#define __problem7_functions_hpp__

#include <iostream>
#include <string>
#include <armadillo>

//Algorithm that solves the thomas algorithm
arma::mat thomas_algo(int n_step, int subdiagonal, int diagonal, int superdiagonal);

// function that writes vectors v and x to file
void write_thomas_to_file(std::string filename, arma::mat v_g_mat);

#endif
