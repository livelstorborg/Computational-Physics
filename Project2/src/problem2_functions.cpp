#include <iostream>
#include <armadillo>
#include <cmath>

arma::mat set_up_A_matrix(const int N){

    /* This function is taking a constant integer N as argument, 
    setting up a tridiagonal matrix A, and returning A*/

    double n = N + 1.0;
    double h = 1.0 / n;
    long double a = - 1.0 / (h * h);
    long double d = 2.0 / (h * h);

    //Setting up A-matrix
    arma::mat A(N, N);

    for (int i = 0; i < N; ++i) {

        A(i, i) = d;            // Maindiagonal

        if (i > 0) {
            A(i, i - 1) = a;    // Subdiagonal
        }

        if (i < N - 1) {
            A(i, i + 1) = a;    // Superdiagonal
        }
    }
    return A;
}

arma::vec solve_eigval(arma::mat A){

    /* This function is taking an armadillo matrix as argument, using arma:eig_sym to compute the eigenvalues 
    and eigenvectors, and returning an armadillo vector with all the eigenvalues*/

    arma::vec eigenvalues;
    arma::mat eigenvectors;

    arma::eig_sym(eigenvalues, eigenvectors, A);

    return eigenvalues;

}

arma::mat solve_eigvec(arma::mat A){

    /* This function is taking an armadillo matrix as argument, using arma:eig_sym to compute the eigenvalues 
    and eigenvectors, and returning an armadillo vector with all the eigenvalues*/

    arma::vec eigenvalues;
    arma::mat eigenvectors;

    arma::eig_sym(eigenvalues, eigenvectors, A);

    return eigenvectors;

}

void analytical_eig_vec_val(arma::mat& analytical_eigvec, arma::vec& analytical_eigval, int N){

    /* This function takes references to an armadillo matrix analytical_eigvec containing the eigen vectors, 
    an armadillo vector and an integer N containing the eigenvalues. 
    This functions the matrix eqation analytically */

    double n = N + 1.0;
    double h = 1.0 / n;
    long double a = - 1.0 / (h * h);
    long double d = 2.0 / (h * h);
    const double pi = M_PI;

    for (int j = 1; j < N + 1; ++j){

        analytical_eigval(j-1) = d + 2 * a * cos(j * pi / (N + 1));     //Analytical lambda

        arma::vec analytical_eigvec_temp(N);                                //Vector for each iteration

        for (int i = 1; i < N + 1; ++i){
            analytical_eigvec_temp(i - 1) = (sin(j * i * pi / (N + 1)));    //Analytical eigenvector
        }
        
        arma::vec analytical_eigvec_temp_norm = arma::normalise(analytical_eigvec_temp);    //Need to normalize the vector
        analytical_eigvec.row(j-1) = analytical_eigvec_temp_norm.t();                       //Transposing and appending
    }
}


void test_eigval_eigvec(int N){

    /* This test function checks that the eigenvalues and eigenvectors from 
    Armadillo agrees with the analytical result for N=6*/

    arma::mat A_analytical = set_up_A_matrix(N);
    arma::mat R = arma::eye(N, N);
    int k = 0;
    int l = 0;

    //Analytical solution
    arma::mat analytical_eigvec(N, N);
    arma::vec analytical_eigval(N);
    analytical_eig_vec_val(analytical_eigvec, analytical_eigval, N);


    //Printing analytical eigenvalues and eigenvectors
    analytical_eigval.print("Eigenvalues analytical:");
    std::cout << "\n";
    analytical_eigvec.print("Eigenvectors analytical:"); //eigenvec 2 and 6 have the 'wrong' sign but that's okay, c_1 = c_4 = -1 is allowed
    std::cout << "\n";

    //Numerical solution
    arma::mat A = set_up_A_matrix(N);
    arma::mat numerical_eigvec = solve_eigvec(A);
    arma::vec numerical_eigval = solve_eigval(A);

    //Printing numerical eigenvalues and eigenvectors
    numerical_eigval.print("Eigenvalues numerical:");
    std::cout << "\n";
    numerical_eigvec.print("Eigenvectors numerical:");
    std::cout << "\n";
    
    //Differences between anaylical and numerical eigvalues and eigvectors
    arma::vec eigval_diff = numerical_eigval - analytical_eigval;                  //Calculates difference in eigenvalues
    arma::mat eigvec_diff(N, N);                                 //Calculates difference in eigenvectors

    //Finds the signs, ref sign differences in eigenvecs 2 and 6
    arma::mat analytical_eigvec_sign = sign(analytical_eigvec); 
    arma::mat numerical_eigvec_sign = sign(numerical_eigvec);

    //Loops over all eigenvectors
    for(int i = 0; i < N; ++i){

        if(analytical_eigvec_sign.col(i)[0] == numerical_eigvec_sign.col(i)[0]){   //Checks if the signs are in agreement
            
            eigvec_diff.col(i) = analytical_eigvec.col(i) - numerical_eigvec.col(i);       //If they are; Subtract
        }

        else{
            
            eigvec_diff.col(i) = analytical_eigvec.col(i) + numerical_eigvec.col(i);       //If they are not; Double negative
        }
    }

    //Printing the difference between numerical and analytical eigenvalues and eigenvectors
    eigval_diff.print("Difference between numerical and analytical eigenvalues:");
    std::cout << "\n";

    eigvec_diff.print("Difference between numerical and analytical eigenvectors:");

}
