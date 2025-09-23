#include <iostream>
#include <armadillo>
#include "problem2_functions.hpp"

int main(int argc, char* argv[]){

    if(argc != 2){
        std::cout << "Provide an integer N for NxN-matrix! "<< std::endl;
        std::cout << "Example: " << argv[0] << " 6" << std::endl;
        return 1;  // Exit with error code
    }

    int N = atoi(argv[1]);


	/* checks that the eigenvalues and eigenvectors from Armadillo agrees 
	with the analytical result for N=6*/
	test_eigval_eigvec(N);

	return 0;
}