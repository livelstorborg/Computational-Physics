#ifndef __problem2_functions_hpp__
#define __problem2_functions_hpp__

#include <iostream>
#include <string>
#include <armadillo>

//function declaration for u(x)
double u(double x);

//function for filling x_vec with equally spaced values
arma::vec fill_x_vec(int start, int stop, int n);

//function that takes x_vec values and creates ux_vec with u(x) values
arma::vec fill_ux_vec(arma::vec x_vec);

//function that writes x and u(x) values to file
void write_to_file(std::string filename, arma::vec x_vec, arma::vec ux_vec);

#endif