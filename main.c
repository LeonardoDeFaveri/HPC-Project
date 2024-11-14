#include<mpi.h>
#include<stdio.h>

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Finalize();
    return 0;
}