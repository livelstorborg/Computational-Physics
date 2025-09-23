# FYS3150 Project1

## File organization:
This project is organized inside the Project1 folder in this manner:
- Project1 folder contains all main.cpp files and all .pdf and .txt will be saved here after running programs.
- src folder: includes all the .cpp and .py scripts
- include folder: includes headers


## How to compile and run this project:

- make (to compile all .cpp files)
- make run_all (to run all executable files)
- make filename.cpp file (to compile specific .cpp file)
- make clean (clean all executables, plots and .txt files)

### Example1: Compiling and running the whole project
- make
- make run_all
- python 3 src/problem2_plotting.py
- python 3 src/problem7_plotting.py
- python 3 src/problem8_plotting.py


### Example2: compiling and running problem 7 with python plotting file
- make problem7_main 
- ./problem7_main 
- python3 src/problem7_plotting.py


# Compiling .cpp files without using make:

## Problem 2
- g++-14 problem2_main.cpp src/problem2_functions.cpp -Iinclude -larmadillo -o problem2_main

## Problem 7
- g++-14 problem7_main.cpp src/problem7_functions.cpp -Iinclude -larmadillo -o problem7_main 

## Problem 8
- g++-14 problem8_main.cpp src/problem8_functions.cpp -Iinclude -larmadillo -o problem8_main

## Problem 9
- g++-14 problem9_main.cpp src/problem9_functions.cpp -Iinclude -larmadillo -o problem9_main

## problem 10
- g++-14 problem10_main.cpp src/problem10_functions.cpp src/problem7_functions.cpp src/problem9_functions.cpp -Iinclude -larmadillo -o problem10_main
