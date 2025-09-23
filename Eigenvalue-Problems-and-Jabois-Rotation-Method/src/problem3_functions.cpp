#include <iostream>
#include <armadillo>
#include <cmath>
#include "problem3_functions.hpp"

//a)
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){ 

    /*This functions takes references to a constant armadillo matrix A, an integer k and an integer l.
    The function identifies the largest off-diagonal element (in absolute values) of the matrix
    under the assumption of a symmetric matrix, and returns an integer offDiag */


    int N = A.n_rows;
    double offDiag = 0.0;

    for (int i = 0; i < N; ++i){

        for (int j = 0; j < N; ++j){

            if (i != j){                                //Ignores the diagonal
                if (fabs(offDiag) < fabs(A(i, j))){     //Checks if the current offDiagonal value of A is larger than a previous one (absolute values)
                    offDiag = A(i, j);                  //New largest offDiagonal value found
                    k = i;
                    l = j;
                }
            }
        }
    }
    return offDiag;
}

//b)
void test_max_offdiag_symmetric(int N){


    /* This function tests the function max_offdiag_symmetric() with a predefined matrix.*/ 

    arma::mat A = arma::eye(N, N);

    A(3,0) = 0.5;
    A(2,1) = -0.7;
    A(1,2) = -0.7;
    A(0,3) = 0.5;

    double offDiag_returned = 0.0; 
    int l;
    int k;

    offDiag_returned = max_offdiag_symmetric(A, l, k); 

    A.print("A-matrix:");
    std::cout << "\n";

    std::cout << "Largest off-diagonal value in A: " << offDiag_returned << " (Absolute values) " << std::endl;

    std::cout << "k: " << k << ", " << " l: " << l << std::endl;
}