# FYS3150 Project2

## File organization:
This project is organized inside the Project2 folder in this manner:
- Project2 folder contains all main.cpp files and all .pdf and .txt will be saved here after running programs.
- src folder: includes all the .cpp and .py scripts
- include folder: includes headers


## How to compile and run this project:

- make (to compile all .cpp files)
- make run_all (to run all executable files)
- make problemx_main (to compile specific .cpp file)
- make clean (clean all executables, plots and .txt files)

- Some files will ask for for an integer N to express an N x N matrix as a command line argument to save time. There will be a warning while running the executable if an N is needed.

### Example1: Compiling and running the whole project
- make
- make run_all (This will run problem 2, 3, 4 with N = 6 and problem 6 with N = 10 and N = 100)
- python3 src/problem6_plotting.py

### Example2: compiling and running problem 7 with python plotting file
- make problem6_main 
- ./problem6_main 
- python3 src/problem6_plotting.py


# Compiling .cpp files without using make:

## problem2_main:
	g++-14 problem2_main.cpp src/problem2_functions.cpp -Iinclude -larmadillo -o problem2_main

## problem3_main:
	g++-14 problem3_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp -Iinclude -larmadillo -o problem3_main

## problem4_main:
	g++-14 problem4_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp src/problem4_functions.cpp -Iinclude -larmadillo -o problem4_main

## problem5_main:
	g++-14 problem5_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp src/problem4_functions.cpp src/problem5_functions.cpp -Iinclude -larmadillo -o problem5_main

## problem6_main:
	g++-14 problem6_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp src/problem4_functions.cpp src/problem5_functions.cpp src/problem6_functions.cpp -Iinclude -larmadillo -o problem6_main
