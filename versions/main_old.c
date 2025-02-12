/**
 * This version of the program models fishes so that each of them has every
 * information required for both fish advancement and experimen setup (e.g.
 * test function, step_ind and step_vol). This version uses one call to
 * `MPI_Allgather` before each operation with which all the information that
 * are necessary for any operation are shared among all fishes. This is
 * inherently ineficient because data not required to compute a movement step
 * or feeding are shared nonetheless.
 * 
 * TESTING:
 * To test this version compiler the program with:
 * `make old`
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "fss_old.h"
#include "../test_functions.h"

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

/**
 * Defines how many times the algorithms is executed for a single configuration
 * (world_size, fishes_count).
 */
#define REP_COUNT 5

void run(
  int world_size, int rank, int total_fishes, struct func_t function,
  MPI_Datatype *mpi_fish_info
);
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
      double elapsed_time = 0;
      for (int j = 0; j < REP_COUNT; j++) {
        // Waits for every process to arrive here before proceeding
        MPI_Barrier(MPI_COMM_WORLD);
        double start_time = MPI_Wtime();
        run(world_size, rank, fishes_count, function, &mpi_fish_info);
        elapsed_time += MPI_Wtime() - start_time;
      }
      // Takes the average of the results
      elapsed_time /= REP_COUNT;
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

void run(
  int world_size, int rank, int total_fishes, struct func_t function,
  MPI_Datatype *mpi_fish_info
) {
  srand(time(NULL) + rank);

  /*
   * Computes the total amount of fishes. The strategy is that of distributing
   * fishes to that all first n-1 processes all have the same number n1 of fishes,
   * process n has n2 < n1 fishes or 0 fishes.
   * This allows for the use of AllGather instead on AllGatherv when the
   * number of fishes cannot be equally distributed among processes. This way,
   * process n still send n1 fishesh, but the one in excess are
   * ghost values not actually part of any computation.
   */
  int rem = total_fishes % world_size;
  int total_local_fishes = total_fishes / world_size;
  int tot = total_local_fishes + 1;
  if (rem != 0) {
    if ((rank + 1) * tot <= total_fishes) {
      total_local_fishes = tot;
    } else {
      // Leave all the remaining fishes to this process. This is guaranteed to
      // be < total_local_fishes of the previous process, because otherwise
      // rem != 0 would have been false
      total_local_fishes = total_fishes - rank * tot;
    }
  } else {
    tot--;
  }

  int total_size = world_size * tot;
  fish_info_t* fishes = malloc(sizeof(fish_info_t) * total_size);
  fish_t* local_fishes = malloc(sizeof(fish_t) * tot);

  // Initialize more-fishes than necessary, but fishes in excess won't be used
  for (int i = 0; i < tot; i++) {
    init(&local_fishes[i], &function);
  }

  // Open file for writing
  FILE *file;
  if (rank == 0) {
    file = fopen("HPC-Project/fish_positions.csv", "w");
    fprintf(file, "cycle,rank,fish_id,position_x,position_y,weight\n");
  }

  /****************************************************************************/
  /***** JUST FOR PLOTTING NECESSITIES, REMOVE FOR PERFORMANCE EVALUATION *****/
  /****************************************************************************/
  #ifdef DEBUG
  for (int i = 0; i < tot; i++) {
    fishes[rank * tot + i] = local_fishes[i].info;
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fishes, tot, *mpi_fish_info, MPI_COMM_WORLD);
  if (rank == 0) {
    int proc = -1;
    for (int i = 0; i < total_fishes; i++) {
      if (i % tot == 0) {
        proc++;
      }
      fprintf(
        file, "%d,%d,%d,%f,%f,%f\n", -1, proc, i,
        fishes[i].positions[0], fishes[i].positions[1], fishes[i].weight
      );
    }
  }
  #endif
  /****************************************************************************/

  for (int cycle = 0; cycle < CYCLES_LIMIT; cycle++) {
    for (int i = 0; i < total_local_fishes; i++) {
      individual_move(&local_fishes[i]);
      PRINT_POS0(
        "After individual move", cycle, rank, i,
        local_fishes[i].info.positions[0], local_fishes[i].info.positions[1],
        local_fishes[i].info.weight
      );
    }

    for (int i = 0; i < tot; i++) {
      fishes[rank * tot + i] = local_fishes[i].info;
    }

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fishes, tot, *mpi_fish_info, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      // This requires `food_improvement` of every fish
      feeding_operator(&local_fishes[i], fishes, total_fishes);
      PRINT_POS0(
        "After feeding operator", cycle, rank, i,
        local_fishes[i].info.positions[0], local_fishes[i].info.positions[1],
        local_fishes[i].info.weight
      );
    }

    for (int i = 0; i < tot; i++) {
      fishes[rank * tot + i] = local_fishes[i].info;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fishes, tot, *mpi_fish_info, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      // This requires `food_improvement` and `displacement` of every fish
      collective_instinctive_move(&local_fishes[i], fishes, total_fishes);
      PRINT_POS0(
        "After collective instinctive move", cycle, rank, i,
        local_fishes[i].info.positions[0], local_fishes[i].info.positions[1],
        local_fishes[i].info.weight
      );
    }

    for (int i = 0; i < tot; i++) {
      fishes[rank * tot + i] = local_fishes[i].info;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fishes, tot, *mpi_fish_info, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      collective_volitive_move(&local_fishes[i], fishes, total_fishes);
      PRINT_POS0(
        "After collective volitive move", cycle, rank, i,
        local_fishes[i].info.positions[0], local_fishes[i].info.positions[1],
        local_fishes[i].info.weight
      );

      decrease_step(&local_fishes[i]);
    }

    /**************************************************************************/
    /**** JUST FOR PLOTTING NECESSITIES, REMOVE FOR PERFORMANCE EVALUATION ****/
    /**************************************************************************/
    #ifdef DEBUG
    for (int i = 0; i < tot; i++) {
      fishes[rank * tot + i] = local_fishes[i].info;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fishes, tot, *mpi_fish_info, MPI_COMM_WORLD);
    if (rank == 0) {
      int proc = -1;
      for (int i = 0; i < total_fishes; i++) {
        if (i % tot == 0) {
          proc++;
        }
        fprintf(
          file, "%d,%d,%d,%f,%f,%f\n", cycle, proc, i,
          fishes[i].positions[0], fishes[i].positions[1], fishes[i].weight
        );
      }
    }
    #endif
    /**************************************************************************/
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