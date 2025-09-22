# FYS3150 Project 5

## File organization
This project is organized inside the Project5 folder in this manner:
- Project5 folder contains a main file needed to run the code for this project, a src folder, an include folder and a folder with csv files; pre_computed.
- src: contains the .cpp file for the class PDEModel, and the .py file needed for simulation and animation in this project. 
- include: contains the headerfile for the .cpp file located in src.
- pre_computed: computationally heavy csv files.

## Compiling and running the whole project (with SuperLu)
If you have the correct verison of SuperLu on your computer, you can compile, run and produce simulations and animations as following 

    make
    make run_all
    python3 src/Plotting.py


## Compiling and running the whole project (without SuperLU)
If the correct version of SuperLu isn´t availible on your system, you can simply use the pre computed csv files stored in the folder pre_computed. 
You then need to modify the Plotting.py file to use the correct path for accessing the pre computed files. 
Modify the first lines in Plotting.py in this manner

    potential_file = 'pre_computed/Potential_2.csv'
    P_file = 'pre_computed/P_2.csv'
    U_real_file = 'pre_computed/U_real_2.csv'
    U_imag_file = 'pre_computed/U_imag_2.csv'
    prob_file = 'pre_computed/Prob_vec_2.csv' 
Followed by 

    python3 src/Plotting.py


### Note 
#### Libraries
To run this code it is necessary to have the Armadillo and SuperLU libraries on your computer.
- Armadillo: for linear algebra operations
- SuperLu: for the direct solution of large, sparse, nonsymmetric systems of linear equations.

The Armadillo solver spsolve can solve matrix equations using either SuperLU or LAPACK. However, in this project, we are dealing with a sparse linear system, where most of the matrix elements are zero. Using LAPACK for such systems is inefficient, as it does not account for sparsity and processes all zero entries, leading to significantly higher computation time and memory usage. Therefore, it is strongly recommended to use SuperLU for solving these equations.

If SuperLU is not available, you can instead use the precomputed solutions provided in the CSV files located in the precomputed folder.

#### Compiler
To compile this project, we used the Clang C++ compiler. If it's not available on your system, you can also compile the project using the GNU C++ compiler (g++).

#### Make clean
To remove all executable files, .csv files, .pdf files and .mp4 files, run the following command

    make clean 
