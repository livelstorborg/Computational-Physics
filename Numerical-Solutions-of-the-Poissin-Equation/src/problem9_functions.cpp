#include "problem9_functions.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <iomanip>

arma::mat special_thomas_algo(int n_step){

	double h = 1.0/n_step;

	arma::vec x = arma::linspace<arma::vec>(0, 1, n_step + 1); //x-vector where x_i is in [0, 1]
	arma::vec f = 100 * exp((-10) * x); //f(x) vector

	arma::vec a = arma::vec(n_step - 2).fill(-1.); //subdiagonal
	arma::vec b = arma::vec(n_step - 1).fill(2.); //diagonal
	arma::vec c = arma::vec(n_step - 2).fill(-2.); //superdiagonal

	arma::vec g = arma::vec(n_step - 1).fill(0.);
	arma::vec v = arma::vec(n_step + 1).fill(0.);

	arma::vec b_tilde = arma::vec(n_step - 1).fill(0);
	arma::vec g_tilde = arma::vec(n_step - 1).fill(0);


	double h_squared = h * h;
	for(int i = 1; i < n_step + 1 ; i++){
		g[i - 1] = h_squared * f[i];
	}

	b_tilde[0] = b[0];
	g_tilde[0] = g[0];

	//Forward substitution -- finding b_tilde and g_tilde values
	//Setting m[i] = a[i] / b[i-1] to save n_step FLOPs
	for(int i = 1; i < n_step - 1; i++){
		double m = (-1) / b_tilde[i - 1];
		b_tilde[i] = 2 - (m * (-2));
		g_tilde[i] = g[i] - m * g_tilde[i-1];
	}

	//backwards substtutuion -- finding
	v[n_step - 1] = g_tilde[n_step - 2] / b_tilde[n_step - 2];
	for(int i = n_step - 2; i >= 1; i--){
		v[i] = (g_tilde[i - 1] - ((-2) * v[i+1]))/b_tilde[i - 1];
	}

	//setting boundary points:
	v[0] = 0;
	v[v.size() - 1] = 0;

	//arranging v and x in matrix:
	arma::mat v_x_mat(n_step + 1, 2);
	v_x_mat.col(0) = v;
	v_x_mat.col(1) = x; 

	return v_x_mat;

}

















