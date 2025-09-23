#include <iostream>
#include <armadillo>
#include <cmath>
#include "problem2_functions.hpp"
#include "problem3_functions.hpp"
#include "problem4_functions.hpp"

// a)
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int& iterations){

    /* This function takes references to an armadillo matrix A and an armadillo matrix R, and two integers as arguments.
    The functions is an implementation of Jacobi´s rotation method.*/

    //Step 1)
    double eps = 1.0e-8;  // Tolerance
    int row = A.n_rows;

    //Step 2) finding k and l 
    double offDiag = max_offdiag_symmetric(A, k, l);
    double a_kl = A(k, l);  // Initial value to start the loop

    //Step 3)
    while(fabs(a_kl) > eps){

        a_kl = A(k, l);
        double a_ll = A(l, l);
        double a_kk = A(k, k);

        //Step 3.1)
        double tau = (a_ll - a_kk) / (2 * a_kl);

        //Step 3.2) Choose the solution that gives the smallest tθ
        double t_theta;

        if(tau > 0){
            t_theta = 1 / (tau + sqrt(1 + tau * tau));
        }

        else{
            t_theta = 1 / (tau - sqrt(1 + tau * tau));
        }

        //Step 3.3
        double c_theta = 1 / sqrt(1 + t_theta * t_theta);
        double s_theta = c_theta * t_theta;

        double a_kk_mp1 = a_kk * (c_theta * c_theta) - (2 * a_kl * c_theta * s_theta) + a_ll * (s_theta * s_theta);
        double a_ll_mp1 = a_ll * (c_theta * c_theta) + (2 * a_kl * c_theta * s_theta) + a_kk * (s_theta * s_theta);

        A(k, k) = a_kk_mp1;
        A(l, l) = a_ll_mp1;
        A(k, l) = 0.0;
        A(l, k) = 0.0;

        // Updating all other elements where i!=k and i!=l
        for(int i = 0; i < row; ++i){

            if(i != k && i != l){

                double a_ik = A(i, k);
                double a_il = A(i, l);

                A(i, k) = a_ik * c_theta - a_il * s_theta;
                A(k, i) = A(i, k);
                A(i, l) = a_il * c_theta + a_ik * s_theta;
                A(l, i) = A(i, l);
            }

            //Step 3.4
            double r_ik = R(i, k);
            double r_il = R(i, l);

            R(i, k) = r_ik * c_theta - r_il * s_theta;
            R(i, l) = r_il * c_theta + r_ik * s_theta;
        }

        //Step 3.5
        offDiag = max_offdiag_symmetric(A, k, l);
        iterations += 1;
    }
    
}

void test_jacobi_rotate(int N){

    /* This function tests the function jacobi_rotate() with a 6x6-matrix.*/

    arma::mat A = set_up_A_matrix(N);
    arma::mat R = arma::eye(N, N);

    int k, l;
    int iterations = 0;

    jacobi_rotate(A, R, k, l, iterations);

    arma::vec numerical_eigval = arma::vec(A.n_rows);
    for(int i = 0; i < A.n_rows; i++){

        numerical_eigval[i] = A(i, i);
    }

    // Analytical solution
    arma::mat analytical_eigvec(N, N);
    arma::vec analytical_eigval(N);
    analytical_eig_vec_val(analytical_eigvec, analytical_eigval, N);

    // sorting values
    arma::uvec analytical_eigval_sorted = arma::sort_index(analytical_eigval);
    analytical_eigval = analytical_eigval(analytical_eigval_sorted);
    analytical_eigvec = analytical_eigvec.cols(analytical_eigval_sorted);

    // Print eigenvalues and eigenvectors
    analytical_eigval.print("Eigenvalues analytical:");
    analytical_eigvec.print("Eigenvectors analytical:");

    // Printing numerical solution
    arma::uvec numerical_eigval_sorted = arma::sort_index(numerical_eigval);
    numerical_eigval = numerical_eigval(numerical_eigval_sorted);
    R = R.cols(numerical_eigval_sorted);

    numerical_eigval.print("Eigenvalues numerical:");
    R.print("Eigenvectors numerical:");

    // Differences
    arma::vec eigval_diff = numerical_eigval - analytical_eigval;
    arma::mat eigvec_diff(N, N);

    arma::mat analytical_eigvec_sign = sign(analytical_eigvec);
    arma::mat numerical_eigvec_sign = sign(R);

    for(int i = 0; i < N; ++i){

        if(analytical_eigvec_sign.col(i)[0] == numerical_eigvec_sign.col(i)[0]){ //Checks if the signs are in agreement
            
            eigvec_diff.col(i) = analytical_eigvec.col(i) - R.col(i);
        }
        else{ //If they are not; Double negative

            eigvec_diff.col(i) = analytical_eigvec.col(i) + R.col(i);
        }
    }

    // Printing the difference between numerical and analytical eigenvalues and eigenvectors
    eigval_diff.print("Difference between numerical and analytical eigenvalues:");
    std::cout << "\n";
    eigvec_diff.print("Difference between numerical and analytical eigenvectors:");
}

