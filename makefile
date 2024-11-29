run:
	gcc main.c -lm -fopenmp && ./a.out

parallel-run:
	gcc new-main.c -lm -fopenmp && ./a.out 2 4 8 16

mpi-run:
	mpicc mpi-main.c -lm && mpirun -np 2 ./a.out
