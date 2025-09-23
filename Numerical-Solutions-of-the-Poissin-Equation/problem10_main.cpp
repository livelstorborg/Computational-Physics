#include <armadillo>
#include <iostream>
#include <iomanip>
#include "problem7_functions.hpp"
#include "problem9_functions.hpp"
#include "problem10_functions.hpp"

int main(){
	arma::ivec n_steps = arma::ivec("10 100 1000 10000 100000 1000000");

	int width = 20;
	int prec = 8;

	std::ofstream ofile;
	ofile.open("time_comparison.txt");

	ofile << std::left << std::setw(width) << "n_steps" << std::setw(width) << "Time general (s)" << std::setw(width) << "Time special (s)" << std::endl; //making header for file

	for(int i = 0; i < n_steps.size(); i++){
		double thomas_general_duration_seconds = timer_general(n_steps[i]);
		double thomas_special_duration_seconds = timer_special(n_steps[i]);

		ofile << std::left << std::setw(width) << std::setprecision(prec) << n_steps[i]
		 	  << std::setw(width) << std::setprecision(prec) << thomas_general_duration_seconds 
		 	  <<std::setw(width) << std::setprecision(prec) << thomas_special_duration_seconds << std::endl;
	}
	return 0;
}
