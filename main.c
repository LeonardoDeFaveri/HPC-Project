#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "fss.h"
#include "test_functions.h"

void run(int world_size, int rank, struct func_t function, MPI_Datatype *mpi_fish_info);
MPI_Datatype register_fish_info_t();

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
    const struct func_t function = get_function((enum func_name)atoi(argv[1]));
    MPI_Datatype mpi_fish_info = register_fish_info_t();

    run(world_size, rank, function, &mpi_fish_info);

    MPI_Type_free(&mpi_fish_info);
  }

  MPI_Finalize();
  return 0;
}

void run(int world_size, int rank, struct func_t function, MPI_Datatype *mpi_fish_info) {
  srand(time(NULL) + rank);
  fish_info_t* fishes = malloc(sizeof(fish_info_t) * world_size);
  fish_t fish;
  init(&fish, &function);

  printf("FISH[%d]: initial value: %f\t initial weight: %f\n", rank, fish.info.value, fish.info.weight);

  for (int cycle = 0; cycle < CYCLES_LIMIT; cycle++) {

    individual_move(&fish);
    
    // Each fish requires:
    // * Maximum value improvement among all fishes to compute its new weight
    // * Value improvement and position displacement of each fish to compute
    //    its new position after collective movement
    // * Position and weight of each fish to compute fish school baricenter
    //    (this value could be computer by just on fish and propagated to the
    //    others)
    fishes[rank] = fish.info;
    MPI_Allgather(&fish.info, 1, *mpi_fish_info, fishes, 1, *mpi_fish_info, MPI_COMM_WORLD);

    feeding_operator(&fish, fishes, world_size);
    collective_instinctive_move(&fish, fishes, world_size);
    collective_volitive_move(&fish, fishes, world_size, rank);

    //// Perché serve?
    ////MPI_Allgather(&fish.info, 1, *mpi_fish_info, fishes, 1, *mpi_fish_info, MPI_COMM_WORLD);
    decrease_step(&fish, cycle);
  }

  printf("FISH[%d]: final value: %f\t final weight: %f\n", rank, fish.info.value, fish.info.weight);
}

MPI_Datatype register_fish_info_t() {
  const int n_fields = 5;
  int block_lengths[] = {1, 1, 1, DIM_COUNT, DIM_COUNT};
  MPI_Datatype field_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Datatype mpi_fish_info;
  MPI_Aint offset[n_fields];

  offset[0] = offsetof(fish_info_t, weight);
  offset[1] = offsetof(fish_info_t, value);
  offset[2] = offsetof(fish_info_t, value_improvement);
  offset[3] = offsetof(fish_info_t, positions);
  offset[4] = offsetof(fish_info_t, displacements);

  MPI_Type_create_struct(n_fields, block_lengths, offset, field_types, &mpi_fish_info);
  MPI_Type_commit(&mpi_fish_info);
  return mpi_fish_info;
}