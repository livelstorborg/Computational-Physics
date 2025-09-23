#include <iostream>
#include "problem6_functions.hpp"

int main(int argc, char* argv[]){

    if(argc != 2){
        std::cout << "Provide an integer N for NxN-matrix! "<< std::endl;
        std::cout << "Example: " << argv[0] << " 6" << std::endl;
        return 1;
    }

    int N = atoi(argv[1]);
    
	write_to_file_numerical(N);

	
	return 0;
}