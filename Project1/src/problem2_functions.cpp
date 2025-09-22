#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <iomanip>
#include "problem2_functions.hpp"

//function for calculating u(x)
double u(double x){

	double func = (1 - ((1 - exp(-10)) * x) - exp(-10 * x));
	return func;
}

//function for filling x_vec with equally spaced values
arma::vec fill_x_vec(int start, int stop, int n){

	arma::vec x_vec = arma::linspace<arma::vec>(start, stop, n);
	return x_vec;
}

//function that takes x_vec values and creates ux_vec with u(x) values
arma::vec fill_ux_vec(arma::vec x_vec){

	arma::vec ux_vec(x_vec.size());
	for(int i = 0; i < x_vec.size(); i++){
	        ux_vec[i] = u(x_vec[i]);
	    }

	return ux_vec;
}

//function that writes x and u(x) values to file
void write_to_file(std::string filename, arma::vec x_vec, arma::vec ux_vec){

    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::scientific << std::setprecision(8);
    ofile << "x values   u(x) values" << "\n";

    for (int i = 0; i < x_vec.size(); i++){
        ofile << x_vec[i] << " " << ux_vec[i] << "\n";
    }

    ofile.close();
}