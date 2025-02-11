/**
 * This version of the program models fishes so that each fish only carries the
 * information it requires to evolve (e.g. weight, position, displacement). Data
 * associated to experiment (e.g. test function, step_ind, step_vol) are put in
 * a different structure that is replicated among all processes instead of all
 * fishes. This version uses call to `MPI_Allgather` to transfer only those data
 * that are stricly necessary to perform a movement step or feeding. There is a
 * call to `MPI_Allgather` for each piece of data, so if an operation requires
 * multiple pieces of information, multiple calls are made.
 * 
 * TESTING:
 * To test this version compiler the program with:
 * `make 3_ag`
 */

#include <float.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "fss_a.h"
#include "../test_functions.h"

#ifdef DEBUG
  #define PRINT(f, ...) printf(f, __VA_ARGS__)
  #define PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight) PRINT("FISH[%d-%d] AT CYCLE[%d] %s POS: %f, %f WEIGHT: %f\n", rank, local_id, cycle, desc, pos_x, pos_y, weight)
  #define PRINT_POS0(desc, cycle, rank, local_id, pos_x, pos_y, weight) if (rank == 0) PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight)
#else
  #define PRINT(f, ...)
  #define PRINT_POS(desc, cycle, rank, local_id, pos_x, pos_y, weight)
  #define PRINT_POS0(desc, cycle, rank, local_id, pos_x, pos_y, weight)
#endif

/**
 * Defines how many times the algorithms is executed for a single configuration
 * (world_size, fishes_count).
 */
#define REP_COUNT 5

void run(
  int world_size, int rank, int total_fishes, struct setup_info_t* setup,
  MPI_Datatype* mpi_dimensions_t
);
/**
 * Returns the maximum value in a vector `values` of length `n`.
 */
double max(double* values, int n);
/**
 * This type allows transferring all dimensions as they were a single value.
 */
MPI_Datatype register_dimensions_t();
/**
 * This function allocates a matrix in which all rows are contiguous in memory.
 * 
 * (This was necessary for some of my previous tries, but maybe can be avoided
 * in this version of the code)
 */
double **allocate_matrix (int rows, int cols);

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Datatype mpi_dimensions_t = register_dimensions_t();

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
      double elapsed_time = 0;
      for (int j = 0; j < REP_COUNT; j++) {
        // Waits for every process to arrive here before proceeding
        MPI_Barrier(MPI_COMM_WORLD);
        double start_time = MPI_Wtime();
        init_setup(&setup, &function);
        run(world_size, rank, fishes_count, &setup, &mpi_dimensions_t);
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
  }

  MPI_Finalize();
  return 0;
}

void run(
  int world_size, int rank, int total_fishes, struct setup_info_t* setup,
  MPI_Datatype* mpi_dimensions_t
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
  fish_t* local_fishes = malloc(sizeof(fish_t) * tot);
  // This array will be reused for both weight and food improvements
  double* improvements = malloc(sizeof(double) * total_size);
  double* weights = malloc(sizeof(double) * total_size);
  // This array will be reused for both positions and displacements
  double** positions = allocate_matrix(total_size, DIM_COUNT);

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
  for (int i = 0; i < tot; i++) {
    int index = rank * tot + i;
    for (int j = 0; j < DIM_COUNT; j++) {
      positions[index][j] = local_fishes[i].positions[j];
    }
    weights[index] = local_fishes[i].weight;
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, *positions, tot, *mpi_dimensions_t, MPI_COMM_WORLD);
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, weights, tot, MPI_DOUBLE, MPI_COMM_WORLD);
  if (rank == 0) {
    int proc = -1;
    for (int i = 0; i < total_fishes; i++) {
      if (i % tot == 0) {
        proc++;
      }
      fprintf(
        file, "%d,%d,%d,%f,%f,%f\n", -1, proc, i,
        positions[i][0], positions[i][1], weights[i]
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

    for (int i = 0; i < tot; i++) {
      improvements[rank * tot + i] = local_fishes[i].food_improvement;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, improvements, tot, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      // This requires `food_improvement` of every fish
      feeding_operator(&local_fishes[i], max(improvements, total_fishes));
      PRINT_POS0(
        "After feeding operator", cycle, rank, i,
        local_fishes[i].positions[0], local_fishes[i].positions[1],
        local_fishes[i].weight
      );
    }

    for (int i = 0; i < tot; i++) {
      int index = rank * tot + i;
      for (int j = 0; j < DIM_COUNT; j++) {
        positions[index][j] = local_fishes[i].displacements[j];
      }
    }
    
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, *positions, tot, *mpi_dimensions_t, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      // This requires `food_improvement` and `displacement` of every fish
      collective_instinctive_move(&local_fishes[i], positions, improvements, total_fishes, setup);
      PRINT_POS0(
        "After collective instinctive move", cycle, rank, i,
        local_fishes[i].positions[0], local_fishes[i].positions[1],
        local_fishes[i].weight
      );
    }

    for (int i = 0; i < tot; i++) {
      int index = rank * tot + i;
      for (int j = 0; j < DIM_COUNT; j++) {
        positions[index][j] = local_fishes[i].positions[j];
      }
      weights[index] = local_fishes[i].weight;
      improvements[index] = local_fishes[i].weight_improvement;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, *positions, tot, *mpi_dimensions_t, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, improvements, tot, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, weights, tot, MPI_DOUBLE, MPI_COMM_WORLD);

    for (int i = 0; i < total_local_fishes; i++) {
      //collective_volitive_move(&local_fishes[i], all_fishes, total_fishes, setup);
      collective_volitive_move(&local_fishes[i], positions, weights, improvements, total_fishes, setup);
      PRINT_POS0(
        "After collective volitive move", cycle, rank, i,
        local_fishes[i].positions[0], local_fishes[i].positions[1],
        local_fishes[i].weight
      );
    }

    decrease_step(setup);

    /**************************************************************************/
    /**** JUST FOR PLOTTING NECESSITIES, REMOVE FOR PERFORMANCE EVALUATION ****/
    /**************************************************************************/
    #ifdef DEBUG
    for (int i = 0; i < tot; i++) {
      int index = rank * tot + i;
      for (int j = 0; j < DIM_COUNT; j++) {
        positions[index][j] = local_fishes[i].positions[j];
      }
      weights[index] = local_fishes[i].weight;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, *positions, tot, *mpi_dimensions_t, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, weights, tot, MPI_DOUBLE, MPI_COMM_WORLD);
    if (rank == 0) {
      int proc = -1;
      for (int i = 0; i < total_fishes; i++) {
        if (i % tot == 0) {
          proc++;
        }
        fprintf(
          file, "%d,%d,%d,%f,%f,%f\n", cycle, proc, i,
          positions[i][0], positions[i][1], weights[i]
        );
      }
    }
    #endif
    /**************************************************************************/
  }

  if (rank == 0) {
    fclose(file);
  }

  free(improvements);
  free(weights);
  // Just one free because the matrix comes from a 1-D array
  free(positions);
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

MPI_Datatype register_dimensions_t() {
  MPI_Datatype mpi_position_t;
  MPI_Type_vector(1, DIM_COUNT, 0, MPI_DOUBLE, &mpi_position_t);
  MPI_Type_commit(&mpi_position_t);
  return mpi_position_t;
}

double **allocate_matrix (int rows, int cols) {
  double  *data   = malloc (rows * cols * sizeof(double));
  double **matrix = malloc (rows * sizeof(double *));
  for (int i = 0; i < rows; i++) {
    matrix[i] = & (data[i * cols]);
  }
  return matrix;
}
