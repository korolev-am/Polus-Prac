Репозиторий для задания по курсу суперкомпьютерное моделирование и технологии

Файл main.c - последовательный код программы. \
Файл omp-main.c - реализация с OpenMP \
Файл mpi-main.c - реализация с MPI \
Файл hybrid-main.c - реализация с MPI + OpenMP \

С помощью файла makefile можно запускать программы локально

Для запусков на polus использовались следущие команды:
Гибридная программа
```
compile:
        module avail && module load SpectrumMPI/10.1.0 && mpicc new-main.c -lm -fopenmp -std=c99
run40x40_1_4:
        bsub -n 1 -oo test40x40_1_4.out -eo test40x40_1_4.err -m "polus-c4-ib" -R "affinity[thread(4,same=core)*1]" OMP_NUM_THREADS=4 mpiexec -n 1 ./a.out
run40x40_2_4:
        bsub -n 2 -oo test40x40_2_4.out -eo test40x40_2_4.err -m "polus-c4-ib" -R "affinity[thread(4,same=core)*2]" OMP_NUM_THREADS=4 mpiexec -n 2 ./a.out
run80x90_2_1:
        bsub -n 2 -W "1:30" -q normal -oo 80x90_2_1.out -eo 80x90_2_1.err -m "polus-c2-ib" -R "affinity[thread(1,same=core)*2]" OMP_NUM_THREADS=1 mpiexec -n 2 ./a.out
run80x90_2_2:
        bsub -n 2 -W "1:30" -q normal -oo 80x90_2_2.out -eo 80x90_2_2.err -m "polus-c4-ib" -R "affinity[thread(2,same=core)*2]" OMP_NUM_THREADS=2 mpiexec -n 2 ./a.out
run80x90_2_4:
        bsub -n 2 -W "1:30" -q normal -oo 80x90_2_4.out -eo 80x90_2_4.err -m "polus-c4-ib" -R "affinity[thread(4,same=core)*2]" OMP_NUM_THREADS=4 mpiexec -n 2 ./a.out
run80x90_2_8:
        bsub -n 2 -W "1:30" -q normal -oo 80x90_2_8.out -eo 80x90_2_8.err -m "polus-c4-ib" -R "affinity[thread(8,same=core)*2]" OMP_NUM_THREADS=8 mpiexec -n 2 ./a.out
run160x180_4_1:
        bsub -n 4 -W "3:00" -q normal -oo 160x180_4_1.out -eo 160x180_4_1.err -m "polus-c4-ib" -R "affinity[thread(1,same=core)*4]" OMP_NUM_THREADS=1 mpiexec -n 4 ./a.out
run160x180_4_2:
        bsub -n 4 -W "3:00" -q normal -oo 160x180_4_2.out -eo 160x180_4_2.err -m "polus-c4-ib" -R "affinity[thread(2,same=core)*4]" OMP_NUM_THREADS=2 mpiexec -n 4 ./a.out
run160x180_4_4:
        bsub -n 4 -W "1:30" -q normal -oo 160x180_4_4.out -eo 160x180_4_4.err -m "polus-c4-ib" -R " affinity[thread(4,same=core)*4]" OMP_NUM_THREADS=4 mpiexec -n 4 ./a.out
run160x180_4_8:
        bsub -n 4 -W "1:30" -q normal -oo 160x180_4_8.out -eo 160x180_4_8.err -m "polus-c4-ib" -R "affinity[thread(8,same=core)*4]" OMP_NUM_THREADS=8 mpiexec -n 4 ./a.out
```
Mpi
```
run1:
        mpicc mpi-main.c -lm && bsub -n 1 -oo test1.out -eo test1.err mpiexec -n 1 ./a.out
run2:
        mpicc mpi-main.c -lm && bsub -n 2 -oo test2.out -eo test2.err -R "affinity[core(2)]" mpiexec -n 2 ./a.out
run4:
        mpicc mpi-main.c -lm && bsub -n 4 -oo test4.out -eo test4.err -R "affinity[core(4)]" mpiexec -n 4 ./a.out
```
OMP:
```
run:
        gcc omp-main.c -lm -fopenmp && bsub -n 1 -oo test1.out -eo test1.err ./a.out 2
```
Последовательный код:
```
run:
        gcc main.c -lm && bsub -oo test1.out -eo test1.err ./a.out
```
