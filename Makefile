SHELL := /bin/bash

main: main.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 ./main.c -o main

clean:
	rm launch.sh.*
