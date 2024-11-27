#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "fss.h"
#include "test_functions.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 2) {
    if (rank == 0) {
      fprintf(stderr, "You should provide a function name as an integer in (0, 4)\n");
    }
  } else {
    struct func_t function = get_function((enum func_name)atoi(argv[1]));
    printf("[%d] choosen function: %d\n", rank, (int)function.name);
  }

  MPI_Finalize();
  return 0;
}