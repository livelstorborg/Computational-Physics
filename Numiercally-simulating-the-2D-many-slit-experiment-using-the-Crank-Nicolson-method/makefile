CC = c++-14
FLAGS = -Iinclude -larmadillo -lsuperlu -lblas -llapack


.PHONY: all clean

all: main

main: main.cpp src/PDEModel.cpp
	$(CC) main.cpp src/PDEModel.cpp $(FLAGS) -o main

run_all: all
	./main &

clean:
	rm -f main
	rm -f Potential_2.csv
	rm -f Prob_vec_2.csv
	rm -f animation.mp4
	rm -f P_2.csv
	rm -f U_real_2.csv
	rm -f U_imag_2.csv
	rm -f Screen_three_slit.pdf
	rm -f time_plot_imag_0_double_slit.pdf
	rm -f time_plot_imag_40_double_slit.pdf
	rm -f time_plot_imag_80_double_slit.pdf
	rm -f time_plot_real_0_double_slit.pdf
	rm -f time_plot_real_40_double_slit.pdf
	rm -f time_plot_real_80_double_slit.pdf
	rm -f time_plot_P_0_double_slit.pdf
	rm -f time_plot_P_40_double_slit.pdf
	rm -f time_plot_P_80_double_slit.pdf
	rm -f two_slits.pdf
	rm -f probabilities_with_pot.pdf
