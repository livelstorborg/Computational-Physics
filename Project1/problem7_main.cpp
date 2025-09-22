#include "problem7_functions.hpp"

int main(){


    //Writing v and x values to file, and plotting exact vs numerical solution for n_steps=10
    arma::mat v_x_mat_10 = thomas_algo(10, -1, 2, -1);
    write_thomas_to_file("problem7_v_x_10steps.txt", v_x_mat_10);

    //Writing v and x values to file, and plotting exact vs numerical solution for n_steps=100
    arma::mat v_x_mat_100 = thomas_algo(100, -1, 2, -1);
    write_thomas_to_file("problem7_v_x_100steps.txt", v_x_mat_100);

    //Writing v and x values to file, and plotting exact vs numerical solution for n_steps=1000
    arma::mat v_x_mat_1000 = thomas_algo(1000, -1, 2, -1);
    write_thomas_to_file("problem7_v_x_1000steps.txt", v_x_mat_1000);

    arma::mat v_x_mat_10000 = thomas_algo(10000, -1, 2, -1);
    write_thomas_to_file("problem7_v_x_10000steps.txt", v_x_mat_10000);

    arma::mat v_x_mat_100000 = thomas_algo(100000, -1, 2, -1);
    write_thomas_to_file("problem7_v_x_100000steps.txt", v_x_mat_100000);

    return 0;
}
