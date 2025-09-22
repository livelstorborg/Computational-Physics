#include "problem2_functions.hpp"

int main(){

    //Reading x and ux values from file and plotting the exact solution of u(x)
    arma::vec x_vec = fill_x_vec(0, 1, 100); // x in [0, 1] 
    arma::vec ux_vec = fill_ux_vec(x_vec);
    write_to_file("x_and_ux_values_output.txt", x_vec, ux_vec);

    return 0;
}