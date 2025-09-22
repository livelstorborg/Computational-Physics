CC = c++-14
FLAGS = -Iinclude -larmadillo
OMPFLAGS = -fopenmp

.PHONY: all clean

all: main main_parallel

main: main.cpp src/IsingModel.cpp
	$(CC) main.cpp src/IsingModel.cpp $(FLAGS) -o main

main_parallel: main_parallel.cpp src/IsingModel_parallel.cpp
	$(CC) main_parallel.cpp src/IsingModel_parallel.cpp $(FLAGS) $(OMPFLAGS) -o main_parallel

run_all: all
	./main &

clean:
	rm -f main
	rm -f main_parallel
	rm -f energy_L20_T1.000000_ordered.txt
	rm -f energy_L20_T1.000000_unordered.txt
	rm -f energy_L20_T2.400000_ordered.txt
	rm -f energy_L20_T2.400000_unordered.txt
	rm -f L2_func_of_temp.txt
	rm -f L20_func_of_temp.txt
	rm -f 2x2_energy_500.pdf
	rm -f 2x2_mag_500.pdf
	rm -f 2x2_Cv_500.pdf
	rm -f 2x2_susc_500.pdf
	rm -f 2x2_energy_50000.pdf
	rm -f 2x2_mag_50000.pdf
	rm -f 2x2_Cv_50000.pdf
	rm -f 2x2_susc_50000.pdf
	rm -f mc_cycles_vs_eps_T10.pdf
	rm -f mc_cycles_vs_eps_T24.pdf
	rm -f eps_distribution_hist.pdf
	rm -f temp_eps.pdf
	rm -f temp_magnetisation.pdf
	rm -f temp_heat_cap.pdf
	rm -f temp_sus.pdf
	rm -f problem9_heat_cap.pdf
	rm -f problem9_susceptibility.pdf 

