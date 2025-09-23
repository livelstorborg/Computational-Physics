#ifndef __PDEModel_hpp__
#define __PDEModel_hpp__

#include <iostream>
#include <armadillo>
#include <vector>
#include <iomanip>


class PDEModel
{
    public:
        double dt;
        double dx;
        double dy;
        double x_c;
        double y_c;
        double sigma_x;
        double sigma_y;
        double p_x;
        double p_y;
        


    // Constructor 
    PDEModel(double dt, double dx, double dy, 
               double x_c, double y_c, 
               double sigma_x, double sigma_y, 
               double p_x, double p_y);

    // Initial conditions u_ij^0
    arma::cx_vec initial_state(int M);

    arma::cx_vec normalised_initial_state(arma::cx_vec u);

    void construct_potential(arma::mat& V, const double v_0, const double M, const int nr_slits, const double thickness, const double centre, const double middle_wall, const double opening);

    int pair_to_single_index(int i, int j, int M);

    std::tuple<arma::sp_cx_mat,arma::sp_cx_mat> construct_A_B(const int M, const double dx, const double dt, const arma::mat V);

    void print_sp_matrix_structure(const arma::sp_cx_mat& A);

    arma::cx_vec crank_nicolson(const arma::sp_cx_mat A, const arma::sp_cx_mat B, arma::cx_vec u);
};

#endif