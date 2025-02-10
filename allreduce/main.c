/**
 * This version of the program puts most of the heavy workload on fishes. If a
 * value need other fishes informations to be computed, these are exchanged and
 * then each fish computes the value.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "fss.h"
#include "../test_functions.h"
#include <float.h>

#ifdef DEBUG
  #define PRINT(f, ...) printf(f, __VA_ARGS__)
  #define PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight) PRINT("FISH[%d-%d] AT CYCLE[%d] %s POS: %f, %f WEIGHT: %f\n", rank, local_id, cycle, desc, pos_x, pos_y, weight)
  #define PRINT_POS0(desc, cycle, rank, local_id, pos_x, pos_y, weight) if (rank == 0) PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight)
#else
  #define PRINT(f, ...)
  #define PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight)
  #define PRINT_POS0(desc, cycle, rank, local_id, pos_x, pos_y, weight)
#endif

void run(int world_size, int rank, struct setup_info_t* setup, int total_fishes);
/**
 * Returns the maximum value in a vector `values` of length `n`.
 */
double max(double* values, int n);

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
    struct setup_info_t setup;
    int max_fishes_count = atoi(argv[2]);

    FILE *output;
    if (rank == 0) {
      output = fopen(argv[3], "a");
    }

    for (int fishes_count = 1; fishes_count <= max_fishes_count; fishes_count *= 2) {
      // Waits for every process to arrive here before proceeding
      MPI_Barrier(MPI_COMM_WORLD);
      double start_time = MPI_Wtime();
      init_setup(&setup, &function);
      run(world_size, rank, &setup, fishes_count);
      double elapsed_time = MPI_Wtime() - start_time;

      if (rank == 0) {
        fprintf(output, "%d,%d,%f\n", world_size, fishes_count, elapsed_time);
      }
    }

    if (rank == 0) {
      fclose(output);
    }
  }

  MPI_Finalize();
  return 0;
}

void run(
  int world_size, int rank, struct setup_info_t* setup, int total_fishes) 
{
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

  fish_t* local_fishes = malloc(sizeof(fish_t) * tot);

  // Initialize more-fishes than necessary, but fishes in excess won't be used
  for (int i = 0; i < tot; i++) {
    init(&local_fishes[i], setup);
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
    if (rank == 0) {
      for (int i = 0; i < total_local_fishes; i++) {
        fprintf(
          file, "%d,%d,%d,%f,%f,%f\n", -1, rank, i,
          local_fishes[i].positions[0], local_fishes[i].positions[1], local_fishes[i].weight
        );
      }
    }
  #endif
  /****************************************************************************/

  for (int cycle = 0; cycle < CYCLES_LIMIT; cycle++) {
    for (int i = 0; i < total_local_fishes; i++) {
      individual_move(&local_fishes[i], setup);
      PRINT_POS0(
        "After individual move", cycle, rank, i,
        local_fishes[i].positions[0], local_fishes[i].positions[1],
        local_fishes[i].weight
      );
    }

    // Compute local maximum food improvement
    double local_max_food_improvement = max(&local_fishes->food_improvement, total_local_fishes);

    // Compute global maximum food improvement with a single allreduce
    double global_max_food_improvement;
    MPI_Allreduce(&local_max_food_improvement, &global_max_food_improvement, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      feeding_operator(&local_fishes[i], global_max_food_improvement);      
      PRINT_POS0(
        "After feeding operator", cycle, rank, i,
        local_fishes[i].positions[0], local_fishes[i].positions[1],
        local_fishes[i].weight
      );
    }
    
    collective_instinctive_move(local_fishes, total_local_fishes, setup);
    collective_volitive_move(local_fishes, total_local_fishes, setup);

    decrease_step(setup);

    /**************************************************************************/
    /**** JUST FOR PLOTTING NECESSITIES, REMOVE FOR PERFORMANCE EVALUATION ****/
    /**************************************************************************/
    #ifdef DEBUG
      if (rank == 0) {
        for (int i = 0; i < total_local_fishes; i++) {
          fprintf(
            file, "%d,%d,%d,%f,%f,%f\n", cycle, rank, i,
            local_fishes[i].positions[0], local_fishes[i].positions[1], local_fishes[i].weight
          );
        }
      }
    #endif
    /**************************************************************************/
  }

  if (rank == 0) {
    fclose(file);
  }

  free(local_fishes);
}

double max(double* values, int n) {
  double m = -DBL_MAX;
  for (int i = 0; i < n; i++) {
    if (values[i] > m) {
      m = values[i];
    }
  }
  return m;
}