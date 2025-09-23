#include <iostream>
#include <armadillo>
#include "problem2_functions.hpp"
#include "problem3_functions.hpp"

int main(int argc, char* argv[]){

    if(argc != 2){
        std::cout << "Provide an integer N for NxN-matrix! "<< std::endl;
        std::cout << "Example: " << argv[0] << " 6" << std::endl;
        return 1;
    }

    int N = atoi(argv[1]);

    test_max_offdiag_symmetric(N);
    
	return 0;
}