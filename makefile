run:
	gcc main.c -lm -fopenmp && ./a.out

omp-run:
	gcc omp-main.c -lm -fopenmp && ./a.out 2

mpi-run:
	mpicc mpi-main.c -lm && mpirun -np 1 ./a.out

hybrid-run:
	mpicc hybrid-main.c -lm -fopenmp && OMP_NUM_THREADS=2 mpiexec -n 2 ./a.out
