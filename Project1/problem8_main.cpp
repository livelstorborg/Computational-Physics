#include "problem8_functions.hpp"
#include "problem7_functions.hpp"
#include <iomanip>
#include <armadillo>



int main(){

int i = 10;
while(i <= pow(10, 7)){

	write_x_v_u(i);
	i *= 10;

	}
return 0;

}