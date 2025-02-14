#include "fss_allreduce.h"
#include <stdbool.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

/******************************************************************************/
/*** UTILITY FUNCTIONS ********************************************************/
/******************************************************************************/
/**
 * Returns a random value between a and b (exclusive).
 */
double rand_real(double min, double max) {
  double div = RAND_MAX / (max - min);
  return min + (rand() / div);
}

void compute_next_position(const fish_t* const fish, const struct func_t* const func, double step_ind, double* next_pos) {
  for (int i = 0; i < DIM_COUNT; i++) {
    double disp = step_ind * rand_real(-1, 1);
    next_pos[i] = fish->positions[i] + disp;

    // A position can't be outside of the search space.
    if (next_pos[i] > func->params.search_space_max) {
      next_pos[i] = func->params.search_space_max;
    } else if (next_pos[i] < func->params.search_space_min) {
      next_pos[i] = func->params.search_space_min;
    }
  }
}

double compute_weight(const fish_t* const fish, double max_food_improvement) {
  // Compute the new weight and make sure its between bounds
  double new_weight = fish->weight;
  if (max_food_improvement != 0) {
    new_weight += fish->food_improvement / max_food_improvement;
    if (new_weight < 1) 
      new_weight = 1;
    else if (new_weight > W_SCALE)
      new_weight = W_SCALE;
  }
  return new_weight;
}

/// The higher the better
double compute_amount_of_food(double value) {
  // The food amount is computed to be inversely proportional to the function's
  // value
  if(value != 0.0)
    return 1.0 / value;
  else
    return DBL_MAX;
}

/**
 * Returns the maximum food_improvement from an array of fishes
 */
double max_food_improvement(fish_t* fishes, int n) {
  double m = -DBL_MAX;
  for (int i = 0; i < n; i++) {
    if (fishes[i].food_improvement > m) {
      m = fishes[i].food_improvement;
    }
  }
  return m;
}

/**
 * Computes the distance between two positions.
 */
double compute_distance(const double* point1, const double* point2) {
  double distance = 0;
  for (int i = 0; i < DIM_COUNT; i++) {
    double diff = point1[i] - point2[i];
    distance += diff * diff;
  }
  return sqrt(distance);
}

/**
 * Gets the best breeding candidate for this process
 */
struct pair_t get_breeding_candidate(
  const fish_t* local_fishes, int local_count, int rank,
  const double* const best_position, int best_index, int worst_index
) {
  double mate_v = -DBL_MAX;
  int mate = -1;
  for (int i = rank * local_count; i < rank * local_count + local_count; i++) {
    if (i != worst_index && i != best_index) {
      int local_i = i % local_count;
      double distance =
        compute_distance(best_position, local_fishes[local_i].positions);
      if (local_fishes[local_i].weight / distance > mate_v) {
        mate_v = local_fishes[local_i].weight / distance;
        mate = i;
      }
    }
  }
  return (struct pair_t) { .index = mate, .value = mate_v };
}

/******************************************************************************/
/*** FISH API & FSS MOVEMENT **************************************************/
/******************************************************************************/

void init_setup(struct setup_info_t* setup, const struct func_t* const f) {
  setup->search_space_width = f->params.search_space_max - f->params.search_space_min;
  setup->step_ind_perc = IND_INIT_PERCENTAGE;
  setup->step_vol_perc = VOL_INIT_PERCENTAGE;
  // NOTE: If here `CYCLES_LIMIT` is replaced with the actual number of algorithm
  // iterations we do, positions change in a much smother way. If the the number
  // of iterations exceedes the cycles count used here, fishes begin to separate
  // on exceeding iterations. Don't know way.
  setup->step_ind_perc_dec = (IND_INIT_PERCENTAGE - IND_FINAL_PERCENTAGE) / CYCLES_LIMIT;
  setup->step_vol_perc_dec = (VOL_INIT_PERCENTAGE - VOL_FINAL_PERCENTAGE) / CYCLES_LIMIT;
  setup->step_ind = setup->search_space_width * setup->step_ind_perc;
  setup->step_vol = setup->search_space_width * setup->step_vol_perc;
  setup->func = *f;
}

void init(fish_t* const fish, const struct setup_info_t* const setup) {
  fish->weight = W_SCALE / 2;
  for (int i = 0; i < DIM_COUNT; i++) {
    fish->positions[i] = rand_real(setup->func.params.init_min, setup->func.params.init_max);
  }
  fish->food_improvement = 0;
  fish->weight_improvement = 0;
}

void individual_move(fish_t* const local_fishes, int local_count, struct setup_info_t* const setup) {
  for (int j = 0; j < local_count; j++) {
    double curr_food = setup->func.f(local_fishes[j].positions, DIM_COUNT);
    double curr_val = compute_amount_of_food(curr_food);

    double next_pos[DIM_COUNT];
    compute_next_position(&local_fishes[j], &setup->func, setup->step_ind, next_pos);

    double next_food = setup->func.f(next_pos, DIM_COUNT);
    double next_val = compute_amount_of_food(next_food);

    local_fishes[j].food_improvement = next_val - curr_val;
    // Checks if the new position is better than the current one
    if (next_val >= curr_val) {
      // The new position is better, so the fish moves
      for (int i = 0; i < DIM_COUNT; i++) {
        local_fishes[j].displacements[i] = next_pos[i] - local_fishes[j].positions[i];
        local_fishes[j].positions[i] = next_pos[i];
      }
      local_fishes[j].value = next_food;
    } else {
      // The new position is worse, so the fish stays in the current position
      for (int i = 0; i < DIM_COUNT; i++) {
        local_fishes[j].displacements[i] = 0;
      }
      local_fishes[j].value = curr_food;
    }
  }
}

void feeding_operator(fish_t* const local_fishes, int local_count) {
  // Compute maximum food improvement among all processes
  double local_max_food_improvement = max_food_improvement(local_fishes, local_count);
  double global_max_food_improvement;
  MPI_Allreduce(&local_max_food_improvement, &global_max_food_improvement, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for (int i = 0; i < local_count; i++) {
    double new_weight = compute_weight(&local_fishes[i], global_max_food_improvement);
    local_fishes[i].weight_improvement = new_weight - local_fishes[i].weight;
    local_fishes[i].weight = new_weight;
  }
}

void collective_instinctive_move(fish_t* const local_fishes, int local_count, struct setup_info_t* const setup) {
  // Pack local contributions: 
  // - first DIM_COUNT slots for local weighted displacements,
  // - last slot for local food improvement
  double local_sum[DIM_COUNT + 1] = {0};
  for (int i = 0; i < local_count; i++) {
    // Use only non-negative food improvements
    double improvement = fmax(local_fishes[i].food_improvement, 0.0);
    for (int j = 0; j < DIM_COUNT; j++) {
      local_sum[j] += local_fishes[i].displacements[j] * improvement;
    }
    local_sum[DIM_COUNT] += improvement;
  }
  
  // Global reduction: sum weighted displacements and improvements across processes
  double global_sum[DIM_COUNT + 1] = {0};
  MPI_Allreduce(local_sum, global_sum, DIM_COUNT + 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  // Only update positions if there was any improvement
  if (global_sum[DIM_COUNT] > 0.0) {
    double move[DIM_COUNT];
    for (int j = 0; j < DIM_COUNT; j++) {
      move[j] = global_sum[j] / global_sum[DIM_COUNT];
    }
    // Apply the computed collective move to every local fish
    for (int i = 0; i < local_count; i++) {
      for (int j = 0; j < DIM_COUNT; j++) {
        local_fishes[i].positions[j] += move[j];
        // Ensure positions remain within the search space bounds
        if (local_fishes[i].positions[j] > setup->func.params.search_space_max)
          local_fishes[i].positions[j] = setup->func.params.search_space_max;
        else if (local_fishes[i].positions[j] < setup->func.params.search_space_min)
          local_fishes[i].positions[j] = setup->func.params.search_space_min;
      }
    }
  }
}

void collective_volitive_move(fish_t* const local_fishes, int local_count, struct setup_info_t* const setup) 
{
  // Pack local contributions:
  // - first DIM_COUNT slots for local weighted displacements,
  // - one slot for local weigth sum,
  // - last slor for local food improvement
  double local_data[DIM_COUNT + 2] = {0};
  for (int i = 0; i < local_count; i++) {
    for (int j = 0; j < DIM_COUNT; j++) {
      local_data[j] += local_fishes[i].positions[j] * local_fishes[i].weight;
    }
    local_data[DIM_COUNT] += local_fishes[i].weight;
    local_data[DIM_COUNT + 1] += local_fishes[i].weight_improvement;
  }

  // Global reduction: sum weighted displacements, weight and weight improvements across processes
  double global_data[DIM_COUNT + 2] = {0};
  MPI_Allreduce(local_data, global_data, DIM_COUNT + 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double global_baricenter[DIM_COUNT];
  double global_total_weight = global_data[DIM_COUNT];
  double global_total_weight_improvement = global_data[DIM_COUNT + 1];

  // Normalize global baricenter
  for (int j = 0; j < DIM_COUNT; j++) {
    global_baricenter[j] = global_data[j] / global_total_weight;
  }
  
  // Determine the movement direction
  int inc = -1;
  if (global_total_weight_improvement < 0) {
    inc = +1;
  }

  // Update each fish's position based on the computed baricenter
  for (int i = 0; i < local_count; i++) {
    double diff[DIM_COUNT];
    double magnitude = 0;
    for (int j = 0; j < DIM_COUNT; j++) {
      diff[j] = local_fishes[i].positions[j] - global_baricenter[j];
      magnitude += diff[j] * diff[j];
    }
    magnitude = sqrt(magnitude);

    // Normalize the difference vector if needed
    if (magnitude > 0) {
      for (int j = 0; j < DIM_COUNT; j++) {
        diff[j] /= magnitude;
      }
    }

    // Update position and ensure it is within bounds
    for (int j = 0; j < DIM_COUNT; j++) {
      local_fishes[i].positions[j] += inc * setup->step_vol * rand_real(0, 1) * diff[j];
      if (local_fishes[i].positions[j] > setup->func.params.search_space_max) {
        local_fishes[i].positions[j] = setup->func.params.search_space_max;
      } else if (local_fishes[i].positions[j] < setup->func.params.search_space_min) {
        local_fishes[i].positions[j] = setup->func.params.search_space_min;
      }
    }
  }
}

// Decrease step_ind and step_vol values
void decrease_step(struct setup_info_t* setup) {
  setup->step_ind_perc -= setup->step_ind_perc_dec;
  setup->step_ind = setup->search_space_width * setup->step_ind_perc;
  setup->step_vol_perc -= setup->step_vol_perc_dec;
  setup->step_vol = setup->search_space_width * setup->step_vol_perc;
}

void breeding(fish_t* const local_fishes, int local_count, int rank) {
  double worst_v = -DBL_MAX, best_v = DBL_MAX;
  int worst, best;
  for (int i = 0; i < local_count; i++) {
    if (local_fishes[i].value < best_v) {
      best_v = local_fishes[i].value;
      best = i;
    }
    if (local_fishes[i].value > worst_v) {
      worst_v = local_fishes[i].value;
      worst = i;
    }
  }

  // Getting the worst of all fishes
  struct pair_t pair_in = { .index = rank * local_count + worst, .value = worst_v };
  struct pair_t worst_p, best_p;
  MPI_Allreduce(&pair_in, &worst_p, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  
  // Getting the best of all fishes
  pair_in.index = rank * local_count + best;
  pair_in.value = best_v;
  MPI_Allreduce(&pair_in, &best_p, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

  double *best_positions;
  bool free_best = false;
  if (best_p.index / local_count == rank) {
    // The process who owns the best fish has to share its position with all
    // the other processes
    best_positions = local_fishes[best_p.index % local_count].positions;
    MPI_Bcast(best_positions, DIM_COUNT, MPI_DOUBLE, rank, MPI_COMM_WORLD);
    pair_in = get_breeding_candidate(
      local_fishes, local_count, rank, best_positions, best_p.index, worst_p.index
    );
  } else {
    // All the other process must get the position of the best fish
    best_positions = malloc(sizeof(double) * DIM_COUNT);
    free_best = true;
    MPI_Bcast(best_positions, DIM_COUNT, MPI_DOUBLE, best_p.index / local_count, MPI_COMM_WORLD);
    pair_in = get_breeding_candidate(
      local_fishes, local_count, rank, best_positions, best_p.index, worst_p.index
    );
  }

  // Let every process know who owns the best candidate
  struct pair_t mate_p;
  MPI_Allreduce(&pair_in, &mate_p, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  int worst_rank = worst_p.index / local_count;
  if (rank == worst_rank) {
    double best_w;
    if (rank != best_p.index / local_count) {
      MPI_Recv(best_positions, DIM_COUNT, MPI_DOUBLE, best_p.index / local_count, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&best_w, 1, MPI_DOUBLE, best_p.index / local_count, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      best_w = local_fishes[best_p.index % local_count].weight;
    }
    double* mate_positions;
    double mate_w;
    bool free_mate = false;
    if (rank != mate_p.index / local_count) {
      mate_positions = malloc(sizeof(double) * DIM_COUNT);
      free_mate = true;
      MPI_Recv(mate_positions, DIM_COUNT, MPI_DOUBLE, mate_p.index / local_count, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&mate_w, 1, MPI_DOUBLE, mate_p.index / local_count, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      mate_positions = local_fishes[mate_p.index % local_count].positions;
      mate_w = local_fishes[mate_p.index % local_count].weight;
    }

    int i = worst_p.index % local_count;
    fish_t* const worst = &local_fishes[i];
    // This process owns the weakest fish
    worst->weight = (best_w + mate_w) / 2;
    for (int j = 0; j < DIM_COUNT; j++) {
      worst->positions[j] = (best_positions[j] + mate_positions[j]) / 2;
    }

    if (free_mate) {
      free(mate_positions);
    }
  } else {
    if (rank == best_p.index / local_count) {
      MPI_Send(best_positions, DIM_COUNT, MPI_DOUBLE, worst_rank, 0, MPI_COMM_WORLD);
      MPI_Send(&local_fishes[best_p.index % local_count].weight, 1, MPI_DOUBLE, worst_rank, 0, MPI_COMM_WORLD);
    }
    if (rank == mate_p.index / local_count) {
      MPI_Send(local_fishes[mate_p.index % local_count].positions, DIM_COUNT, MPI_DOUBLE, worst_rank, 0, MPI_COMM_WORLD);
      MPI_Send(&local_fishes[mate_p.index % local_count].weight, 1, MPI_DOUBLE, worst_rank, 0, MPI_COMM_WORLD);
    }
  }
  if (free_best) {
    free(best_positions);
  }
}

/******************************************************************************/
