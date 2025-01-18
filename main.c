#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "fss.h"
#include "test_functions.h"

#ifdef DEBUG
  #define PRINT(f, ...) printf(f, __VA_ARGS__)
  #define PRINT_INFO(desc, rank, local_id, info) PRINT("FISH[%d-%d] %s: {\n\tfood: %f\n\tweight: %f\n\tpositions: %f, %f\n}\n", rank, local_id, desc, info.food_amount, info.weight, info.positions[0], info.positions[1])
  #define PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight) PRINT("FISH[%d-%d] AT CYCLE[%d] %s POS: %f, %f WEIGHT: %f\n", rank, local_id, cycle, desc, pos_x, pos_y, weight)
  #define PRINT_POS0(desc, cycle, rank, local_id, pos_x, pos_y, weight) if (rank == 0) PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight)
#else
  #define PRINT(f, ...)
  #define PRINT_INFO(desc, rank, local_id, info)
  #define PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight)
  #define PRINT_POS0(desc, cycle, rank, local_id, pos_x, pos_y, weight)
#endif

void run(int world_size, int rank, struct func_t function, MPI_Datatype *mpi_fish_info, int total_fishes);
MPI_Datatype register_fish_info_t();

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 4) {
    if (rank == 0) {
      fprintf(stderr, "You should provide a function name as an integer in (0, 4), the maximum number of fishes and a path to the output file\n");
    }
  } else {
    const struct func_t function = get_function((enum func_name) atoi(argv[1]));
    int max_fishes_count = atoi(argv[2]);
    MPI_Datatype mpi_fish_info = register_fish_info_t();

    FILE *output;
    if (rank == 0) {
      output = fopen(argv[3], "a");
    }
    for (int fishes_count = 1; fishes_count <= max_fishes_count; fishes_count *= 2) {
      // Waits for every process to arrive here before proceeding
      MPI_Barrier(MPI_COMM_WORLD);
      double start_time = MPI_Wtime();
      run(world_size, rank, function, &mpi_fish_info, max_fishes_count);
      double elapsed_time = MPI_Wtime() - start_time;

      if (rank == 0) {
        fprintf(output, "%d,%d,%f\n", world_size, fishes_count, elapsed_time);
      }
    }
    if (rank == 0) {
      fclose(output);
    }

    MPI_Type_free(&mpi_fish_info);
  }

  MPI_Finalize();
  return 0;
}

void run(int world_size, int rank, struct func_t function, MPI_Datatype *mpi_fish_info, int total_fishes) {
  srand(time(NULL) + rank);

  // Compute the amount of local fishes
  int total_local_fishes = total_fishes / world_size;
  if (rank < total_fishes % world_size) {
    total_local_fishes++;
  }

  fish_info_t* fishes = malloc(sizeof(fish_info_t) * total_fishes);
  fish_info_t* send_buffer = malloc(sizeof(fish_info_t) * total_local_fishes);
  fish_t* local_fishes = malloc(sizeof(fish_t) * total_local_fishes);

  for (int i = 0; i < total_local_fishes; i++) {
    init(&local_fishes[i], &function);
  }

  // Open file for writing
  FILE *file;
  if (rank == 0) {
    file = fopen("HPC-Project/fish_positions.csv", "w");
    fprintf(file, "cycle,rank,fish_id,position_x,position_y,weight\n");
  }

  for (int cycle = 0; cycle < 100; cycle++) {
    for (int i = 0; i < total_local_fishes; i++) {
      individual_move(&local_fishes[i]);
      PRINT_POS0("After individual move", cycle, rank, i, local_fishes[i].info.positions[0], local_fishes[i].info.positions[1], local_fishes[i].info.weight);
    }

    for (int i = 0; i < total_local_fishes; i++) {
      fishes[rank * total_local_fishes + i] = local_fishes[i].info;
    }

    for (int i = 0; i < total_local_fishes; i++) {
      send_buffer[i] = local_fishes[i].info;
    }

    MPI_Allgather(send_buffer, total_local_fishes, *mpi_fish_info, fishes, total_local_fishes, *mpi_fish_info, MPI_COMM_WORLD);

    if (rank == 0) {
      int fish_index = 0;
      for (int proc = 0; proc < world_size; proc++) {
        int fishes_in_proc = total_fishes / world_size;
        if (proc < total_fishes % world_size) {
          fishes_in_proc++;
        }
        for (int i = 0; i < fishes_in_proc; i++) {
          fprintf(file, "%d,%d,%d,%f,%f,%f\n", cycle, proc, fish_index, fishes[fish_index].positions[0], fishes[fish_index].positions[1], fishes[fish_index].weight);
          fish_index++;
        }
      }
    }

    for (int i = 0; i < total_local_fishes; i++) {
      feeding_operator(&local_fishes[i], fishes, total_fishes);
      PRINT_POS0("After feeding operator", cycle, rank, i, local_fishes[i].info.positions[0], local_fishes[i].info.positions[1], local_fishes[i].info.weight);

      collective_instinctive_move(&local_fishes[i], fishes, total_fishes);
      PRINT_POS0("After collective instinctive move", cycle, rank, i, local_fishes[i].info.positions[0], local_fishes[i].info.positions[1], local_fishes[i].info.weight);

      collective_volitive_move(&local_fishes[i], fishes, total_fishes, rank);
      PRINT_POS0("After collective volitive move", cycle, rank, i, local_fishes[i].info.positions[0], local_fishes[i].info.positions[1], local_fishes[i].info.weight);

      decrease_step(&local_fishes[i], cycle);
    }
  }

  if (rank == 0) {
    fclose(file);
  }

  free(fishes);
  free(local_fishes);
}

MPI_Datatype register_fish_info_t() {
  const int n_fields = 6;
  int block_lengths[] = {1, 1, 1, 1, DIM_COUNT, DIM_COUNT};
  MPI_Datatype field_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Datatype mpi_fish_info;
  MPI_Aint offset[n_fields];

  offset[0] = offsetof(fish_info_t, weight);
  offset[1] = offsetof(fish_info_t, food_amount);
  offset[2] = offsetof(fish_info_t, food_improvement);
  offset[3] = offsetof(fish_info_t, weight_improvement);
  offset[4] = offsetof(fish_info_t, positions);
  offset[5] = offsetof(fish_info_t, displacements);

  MPI_Type_create_struct(n_fields, block_lengths, offset, field_types, &mpi_fish_info);
  MPI_Type_commit(&mpi_fish_info);
  return mpi_fish_info;
}