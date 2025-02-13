/**
 * This version of the program models fishes so that each fish only carries the
 * information it requires to evolve (e.g. weight, position, displacement). Data
 * associated to experiment (e.g. test function, step_ind, step_vol) are put in
 * a different structure that is replicated among all processes instead of all
 * fishes. This version uses call to `MPI_Allgather` to transfer only those data
 * that are stricly necessary to perform a movement step or feeding. There is
 * only one call to `MPI_Allgather` for each operation. MPI primitives and
 * derived data types are used to achieve this.
 * 
 * TESTING:
 * To test this version compiler the program with:
 * `make 1_ag`
 */

#include <math.h>
#include <float.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "fss_b.h"
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
  int world_size, int rank, int total_fishes, double* baricenter,
  struct setup_info_t* setup, MPI_Datatype* mpi_dimensions_t,
  MPI_Datatype* mpi_volitive_t
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
 * Packs together all the info necessary to perform the volitive step of a fish.
 */
MPI_Datatype register_volitive_t();
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
  MPI_Datatype mpi_volitive_t = register_volitive_t();

  if (argc < 4) {
    if (rank == 0) {
      fprintf(stderr, "You should provide a function name as an integer in (0, 4), the maximum number of fishes and a path to the output file\n");
    }
  } else {
    const struct func_t function = get_function((enum func_name) atoi(argv[1]));
    struct setup_info_t setup;
    int max_fishes_count = atoi(argv[2]);
    double *baricenter = malloc(sizeof(double) * DIM_COUNT);

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
        run(world_size, rank, fishes_count, baricenter, &setup, &mpi_dimensions_t, &mpi_volitive_t);
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
    free(baricenter);
  }

  MPI_Type_free(&mpi_dimensions_t);
  MPI_Type_free(&mpi_volitive_t);
  MPI_Finalize();
  return 0;
}

void run(
  int world_size, int rank, int total_fishes, double* baricenter,
  struct setup_info_t* setup, MPI_Datatype* mpi_dimensions_t,
  MPI_Datatype* mpi_volitive_t
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
  // Number of real fishes allocated to this process
  int total_local_fishes = total_fishes / world_size;
  // Number of fishes (both real and ghosts) allocated to this process
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
  fish_t* all_fishes = malloc(sizeof(fish_t) * total_size);
  double* food_improvements = malloc(sizeof(double) * total_size);
  double** displacements = allocate_matrix(total_size, DIM_COUNT);

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
    all_fishes[rank * tot + i] = local_fishes[i];
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, all_fishes, tot, *mpi_volitive_t, MPI_COMM_WORLD);
  if (rank == 0) {
    int proc = -1;
    for (int i = 0; i < total_fishes; i++) {
      if (i % tot == 0) {
        proc++;
      }
      fprintf(
        file, "%d,%d,%d,%f,%f,%f\n", -1, proc, i,
        all_fishes[i].positions[0], all_fishes[i].positions[1], all_fishes[i].weight
      );
    }
  }
  #endif
  /****************************************************************************/

  for (int cycle = 0; cycle < CYCLES_LIMIT; cycle++) {
    for (int i = 0; i < total_local_fishes; i++) {
      individual_move(&local_fishes[i], setup);
    }

    for (int i = 0; i < tot; i++) {
      food_improvements[rank * tot + i] = local_fishes[i].food_improvement;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, food_improvements, tot, MPI_DOUBLE, MPI_COMM_WORLD);    
    double max_f = max(food_improvements, total_fishes);
    for (int i = 0; i < total_local_fishes; i++) {
      feeding_operator(&local_fishes[i], max_f);
    }

    for (int i = 0; i < tot; i++) {
      int index = rank * tot + i;
      for (int j = 0; j < DIM_COUNT; j++) {
        displacements[index][j] = local_fishes[i].displacements[j];
      }
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, *displacements, tot, *mpi_dimensions_t, MPI_COMM_WORLD);
    for (int i = 0; i < total_local_fishes; i++) {
      collective_instinctive_move(&local_fishes[i], displacements, food_improvements, total_fishes, setup);
    }

    MPI_Allgather(local_fishes, tot, *mpi_volitive_t, all_fishes, tot, *mpi_volitive_t, MPI_COMM_WORLD);
    compute_baricenter(baricenter, all_fishes, total_fishes);
    double total_weight_improvement = 0;
    for (int i = 0; i < total_fishes; i++) {
      total_weight_improvement += all_fishes[i].weight_improvement;
    }
    for (int i = 0; i < total_local_fishes; i++) {
      collective_volitive_move(&local_fishes[i], baricenter, total_weight_improvement, setup);
    }

    decrease_step(setup);

    // Breeding operator
    double min_v = DBL_MAX, max_v1 = -DBL_MAX, max_v2 = -DBL_MAX;
    int min, max_1, max_2;
    // Gathers the value of the fitness function of each fish
    for (int i = 0; i < tot; i++) {
      food_improvements[rank * tot + i] = local_fishes[i].value;
    }
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, food_improvements, tot, MPI_DOUBLE, MPI_COMM_WORLD);    
    for (int i = 0; i < total_fishes; i++) {
      double value = food_improvements[i];
      if (value > max_v1) {
        max_v2 = max_v1;
        max_2 = max_1;
        max_v1 = value;
        max_1 = i;
      } else if (value > max_v2) {
        max_v2 = value;
        max_2 = i;
      } else if (value < min_v) {
        min_v = value;
        min = i;
      }
    }
    if (rank * tot <= min && (rank + 1) * tot > min) {
      // This process owns the weakest fish
      local_fishes[min].weight = (all_fishes[max_1].weight + all_fishes[max_2].weight) / 2;
      for (int j = 0; j < DIM_COUNT; j++) {
        local_fishes[min].positions[j] = (
          all_fishes[max_1].positions[j] +
          all_fishes[max_2].positions[j]
        ) / 2;
      }
    }

    /**************************************************************************/
    /**** JUST FOR PLOTTING NECESSITIES, REMOVE FOR PERFORMANCE EVALUATION ****/
    /**************************************************************************/
    #ifdef DEBUG
    MPI_Allgather(local_fishes, tot, *mpi_volitive_t, all_fishes, tot, *mpi_volitive_t, MPI_COMM_WORLD);
    if (rank == 0) {
      int proc = -1;
      for (int i = 0; i < total_fishes; i++) {
        if (i % tot == 0) {
          proc++;
        }
        fprintf(
          file, "%d,%d,%d,%f,%f,%f\n", cycle, proc, i,
          all_fishes[i].positions[0], all_fishes[i].positions[1], all_fishes[i].weight
        );
      }
    }
    #endif
    /**************************************************************************/
  }

  if (rank == 0) {
    fclose(file);
  }

  free(food_improvements);
  // Just one free because the matrix comes from a 1-D array
  free(displacements);
  free(local_fishes);
  free(all_fishes);
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

MPI_Datatype register_volitive_t() {
  const int n_fields = 3;
  int block_lengths[] = {DIM_COUNT, 1, 1};
  MPI_Datatype field_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offset[n_fields];
  offset[0] = offsetof(fish_t, positions);
  offset[1] = offsetof(fish_t, weight);
  offset[2] = offsetof(fish_t, weight_improvement);

  MPI_Datatype tmp_mpi_t;
  MPI_Type_create_struct(n_fields, block_lengths, offset, field_types, &tmp_mpi_t);

  /**
   * The following three lines are necessary if we want to send this data type
   * multiple times with the same call (we want this). What these lines do is
   * taking padding and any spaces the the compiler puts in the memory
   * representation of the data type. In this case, since we're sending only
   * some of the struct fields we need to consider that there are other fields
   * in memory.
   */
  MPI_Datatype mpi_volitive_t;
  MPI_Type_create_resized(tmp_mpi_t, offset[0], (MPI_Aint) sizeof(fish_t), &mpi_volitive_t);
  MPI_Type_commit(&mpi_volitive_t);
  return mpi_volitive_t;
}

double **allocate_matrix(int rows, int cols) {
  double  *data   = malloc (rows * cols * sizeof(double));
  double **matrix = malloc (rows * sizeof(double *));
  for (int i = 0; i < rows; i++) {
    matrix[i] = &(data[i * cols]);
  }
  return matrix;
}
