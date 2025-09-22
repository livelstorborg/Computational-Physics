# FYS3150 Project3

## File organization:
This project is organized inside the Project3 folder in this manner:
- Project3 folder contains all main.cpp files and all computationally heavy .txt files will be saved in the folder txt_files after running the programs.
- Other non-computationally heavy txt files will be saved in the Project 3 folder as well.
- src folder: includes all the .cpp and .py scripts
- include folder: includes headers


## How to compile and run this project:

- make (to compile all main files)
- make run_all (to run all executable files)
    - including main.cpp and main_time_dependent.cpp
- make example.cpp (to compile specific .cpp file)
- To make main_time_dependent, you first have to run ./main because main_time_dependent is also dependent on the .txt files made my main.cpp
- make clean (clean all executables, plots and .txt files except the computationally heavy files located in the txt_files folder.)

### Example1: Compiling and running the whole project (however it is strongly recommended to use provided .txt files for plotting instead of making your own by running main_time_dependent.cpp!)

- make
- make run_all
- python3 src/plotting.py

### Example2: compiling and running individually

- make main
- ./main
- python3 src/plotting.py


# Compiling .cpp files without using make:

## main:
    g++-14 main.cpp src/Particle.cpp src/PenningTrap.cpp -Iinclude -larmadillo -o main

## main_time_dependent:
    g++-14 main_time_dependent.cpp src/Particle.cpp src/PenningTrap.cpp -Iinclude -larmadillo -o main -O3
