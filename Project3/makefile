CC = g++-14
IFLAG = -Iinclude

.PHONY: all clean

all: main main_time_dependent

main: main.cpp src/Particle.cpp src/PenningTrap.cpp
	$(CC) main.cpp src/Particle.cpp src/PenningTrap.cpp $(IFLAG) -larmadillo -o main

main_time_dependent: main_time_dependent.cpp src/Particle.cpp src/PenningTrap.cpp
	$(CC) main_time_dependent.cpp src/Particle.cpp src/PenningTrap.cpp $(IFLAG) -larmadillo -o main_time_dependent -O3

run_all: all
	./main &
	./main_time_dependent &


clean:
	rm -f main
	rm -f main_time_dependent
	rm -f xyz_analytical.txt
	rm -f xyz_RK4_1.txt
	rm -f xyz_RK4_1_velocity.txt
	rm -f xyz_RK4_2.txt
	rm -f xyz_RK4_2_velocity.txt
	rm -f xyz_RK4_interactions_1.txt
	rm -f xyz_RK4_interactions_1_velocity.txt
	rm -f xyz_RK4_interactions_2.txt
	rm -f xyz_RK4_interactions_2_velocity.txt
	rm -f xyz_FE_1.txt
	rm -f xyz_FE_1_interactions.txt
	rm -f xyz_FE_2.txt
	rm -f xyz_FE_2_interactions.txt
	rm -f relative_error_FE_4000.txt
	rm -f relative_error_FE_8000.txt
	rm -f relative_error_FE_16000.txt
	rm -f relative_error_FE_32000.txt
	rm -f relative_error_RK4_4000.txt
	rm -f relative_error_RK4_8000.txt
	rm -f relative_error_RK4_16000.txt
	rm -f relative_error_RK4_32000.txt
	rm -f 8.1_z_t.pdf
	rm -f XY-Plane_without_interactions.pdf
	rm -f XY-Plane_analytical.pdf
	rm -f 3D.pdf
	rm -f 3D_interactions.pdf
	rm -f log_error.pdf
	rm -f x_vx.pdf
	rm -f x_vx_interactions.pdf
	rm -f z_vz.pdf
	rm -f z_vz_interactions.pdf
	rm -f XY-Plane_FE_interactions.pdf
	rm -f XY-Plane_RK4_interactions.pdf
	rm -f difference_error_plot.pdf
	rm -f Particle_in_trap_at_omega_v.pdf
	rm -f relative_error.pdf
	rm -f Particle_f_01.pdf
	rm -f Particle_f_04.pdf
	rm -f Particle_f_07.pdf
	rm -f Particle_f_07_both.pdf
	rm -f Particle_f_01_both.pdf

