#include <iostream>
#include "problem2_functions.hpp"
#include "problem3_functions.hpp"
#include "problem4_functions.hpp"

int main(int argc, char* argv[]){

    if(argc != 2){
        std::cout << "Provide an integer N for NxN-matrix! "<< std::endl;
        std::cout << "Example: " << argv[0] << " 6" << std::endl;
        return 1;
    }

    int N = atoi(argv[1]);

	/*b) Testing jacobi_rotate() with N=6 and comparing numerical and analytical results
	by printing to terminal.*/

	test_jacobi_rotate(N);

	return 0;
}