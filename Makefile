SHELL := /bin/bash
#-------------------------------------------------------------------------------
# Comment this line before computing execution times
DEBUG := "-DDEBUG"
#-------------------------------------------------------------------------------

main: main.c test_functions.h test_functions.c fss.h fss.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) ./main.c ./test_functions.c ./fss.c -o main

old: main_old.c test_functions.h test_functions.c fss_old.h fss_old.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) ./main_old.c ./test_functions.c ./fss_old.c -o main

allreduce: allreduce/main.c test_functions.h test_functions.c allreduce/fss.h allreduce/fss.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) ./allreduce/main.c ./test_functions.c ./allreduce/fss.c -o main

.PHONY: clean
clean:
	@rm launch.sh.*
