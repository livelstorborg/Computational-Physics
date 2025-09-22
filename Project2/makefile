CC = g++-14
IFLAG = -Iinclude

.PHONY: all clean

all: problem2_main problem3_main problem4_main problem5_main problem6_main

problem2_main:
	$(CC) problem2_main.cpp src/problem2_functions.cpp $(IFLAG) -larmadillo -o problem2_main

problem3_main:
	$(CC) problem3_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp $(IFLAG) -larmadillo -o problem3_main

problem4_main:
	$(CC) problem4_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp src/problem4_functions.cpp $(IFLAG) -larmadillo -o problem4_main

problem5_main:
	$(CC) problem5_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp src/problem4_functions.cpp src/problem5_functions.cpp $(IFLAG) -larmadillo -o problem5_main

problem6_main:
	$(CC) problem6_main.cpp src/problem2_functions.cpp src/problem3_functions.cpp src/problem4_functions.cpp src/problem5_functions.cpp src/problem6_functions.cpp $(IFLAG) -larmadillo -o problem6_main

run_all: all
	./problem2_main 6 &
	./problem3_main 6 &
	./problem4_main 6 &
	./problem5_main &
	./problem6_main 10 &
	./problem6_main 100 &

clean:
	rm -f problem2_main
	rm -f problem3_main
	rm -f problem4_main
	rm -f problem5_main
	rm -f problem6_main
	rm -f N_vs_similarity_transformation.txt
	rm -f numerical_smallest_lambda.txt
	rm -f numerical_vs_analytical_10.pdf
	rm -f numerical_vs_analytical_100.pdf
