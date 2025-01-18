#!/bin/bash

#PBS -l select=1:ncpus=4:mem=2gb
#PBS -l walltime=0:10:00
#PBS -q short_cpuQ
module load mpich-3.2

# Resets results file
echo "WORLD_SIZE,FISHES_COUNT,ELAPSED_TIME" > HPC-Project/time_results.csv
test_function=0
max_fishes_count=64

# Arguments:
# $1: test function selector
# $2: maximum number of fishes
# $3: files for time_results
mpirun.actual -n 1 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv
mpirun.actual -n 2 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv
mpirun.actual -n 4 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv
mpirun.actual -n 8 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv
mpirun.actual -n 16 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv
mpirun.actual -n 32 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv
mpirun.actual -n 64 ./HPC-Project/main $test_function $max_fishes_count HPC-Project/time_results.csv