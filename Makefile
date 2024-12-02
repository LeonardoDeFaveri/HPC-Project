SHELL := /bin/bash

main: main.c test_functions.h test_functions.c fss.h fss.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl ./main.c ./test_functions.c ./fss.c -o main

clean:
	rm launch.sh.*
