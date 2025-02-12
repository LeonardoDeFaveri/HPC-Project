SHELL := /bin/bash
#-------------------------------------------------------------------------------
# Comment this line before computing execution times
DEBUG := "-DDEBUG"
#-------------------------------------------------------------------------------

main: main.c test_functions.h test_functions.c fss.h fss.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) ./main.c ./test_functions.c ./fss.c -o main

old: versions/main_old.c test_functions.h test_functions.c versions/fss_old.h versions/fss_old.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) versions/main_old.c ./test_functions.c versions/fss_old.c -o main

3_ag: versions/main_a.c test_functions.h test_functions.c versions/fss_a.h versions/fss_a.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) versions/main_a.c ./test_functions.c versions/fss_a.c -o main

allreduce: versions/main_allreduce.c test_functions.h test_functions.c versions/fss_allreduce.h versions/fss_allreduce.c
	@module load mpich-3.2 && mpicc -Wall -std=c11 -lm -ldl $(DEBUG) ./versions/main_allreduce.c ./test_functions.c ./versions/fss_allreduce.c -o main

.PHONY: clean
clean:
	@rm launch.sh.*
