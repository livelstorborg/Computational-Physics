CC = g++-14 
IFLAG = -Iinclude

.PHONY: all clean 

all: problem2_main problem7_main problem8_main problem9_main problem10_main

problem2_main: 
	$(CC) problem2_main.cpp src/problem2_functions.cpp $(IFLAG) -larmadillo -o problem2_main

problem7_main: 
	$(CC) problem7_main.cpp src/problem7_functions.cpp $(IFLAG) -larmadillo -o problem7_main

problem8_main: 
	$(CC) problem8_main.cpp src/problem8_functions.cpp src/problem7_functions.cpp $(IFLAG) -larmadillo -o problem8_main

problem9_main: 
	$(CC) problem9_main.cpp src/problem9_functions.cpp $(IFLAG) -larmadillo -o problem9_main

problem10_main: 
	$(CC) problem10_main.cpp src/problem10_functions.cpp src/problem7_functions.cpp src/problem9_functions.cpp $(IFLAG) -larmadillo -o problem10_main

run_all: all
	./problem2_main &
	./problem7_main &
	./problem8_main &
	./problem9_main	&
	./problem10_main &

clean:
	rm -f problem2_main
	rm -f problem7_main
	rm -f problem8_main
	rm -f problem9_main
	rm -f problem10_main
	rm -f problem7_v_x_10steps.txt
	rm -f problem7_v_x_100steps.txt
	rm -f problem7_v_x_1000steps.txt
	rm -f problem7_v_x_10000steps.txt
	rm -f problem7_v_x_100000steps.txt
	rm -f x_v_u10.txt
	rm -f x_v_u100.txt
	rm -f x_v_u1000.txt
	rm -f x_v_u10000.txt
	rm -f x_v_u100000.txt
	rm -f x_v_u1000000.txt
	rm -f x_v_u10000000.txt
	rm -f log10_relative_error_h.pdf
	rm -f log10_absolute_error_h.pdf
	rm -f x_and_ux_values_output.txt
	rm -f plot_comparing_exact_numerical.pdf
	rm -f problem2_ux_plot.pdf
	rm -f x_log_error_plot.pdf
	rm -f log10_absolute_relative_error_h.pdf
	rm -f max_rel_err_table.txt
	rm -f time_comparison.txt
