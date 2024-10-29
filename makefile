run:
	gcc main.c -lm -fopenmp && ./a.out

parallel-run:
	gcc new-main.c -lm -fopenmp && ./a.out 2 4 8 16