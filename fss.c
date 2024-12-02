#include "fss.h"

#include <time.h>
#include <stdlib.h>

/******************************************************************************/
/*** UTILITY FUNCTIONS *********************************************************/
/******************************************************************************/
/**
 * Returns a random `double` value between a and b (inclusive).
 */
double rand_real(double min, double max) {
  double range = (max - min) + 1; 
  double div = RAND_MAX / range;
  return min + (rand() / div);
}
/**
 * Randomly generates a direction for the displacement over one dimension.
 */
int8_t direction() {
  int8_t v = rand() % 2;
  if (v == 0) {
    return -1;
  }
  return 1;
}

void compute_next_position(const fish_t* fish, double next_pos[]) {
  for (int i = 0; i < DIM_COUNT; i++) {
    double disp = fish->step_ind * rand_real(0, 1) * direction();
    next_pos[i] = fish->info.positions[i] + disp;

    // A position can't be outside of the search space.
    if (
      next_pos[i] > fish->func.params.search_space_max ||
      next_pos[i] < fish->func.params.search_space_min
    ) {
      next_pos[i] = fish->info.positions[i];
    }
  }
}

double compute_weight(const fish_t* fish, fish_info_t* fishes, int n) {
  double max = 0;
  for (int i = 0; i < n; i++) {
    if (fishes[i].value_improvement > max) {
      max = fishes[i].value_improvement;
    }
  }
  return fish->info.weight + fish->info.value_improvement / max;
}

double decrease_linearly(double value, double initial_value, double final_value, int max_iterations) {
  return value - (initial_value - final_value) / max_iterations;
}

/******************************************************************************/
/*** FISH API & FSS MOVEMENT **************************************************/
/******************************************************************************/

void init(fish_t* fish, const struct func_t* func) {
  double width = func->params.search_space_max - func->params.search_space_min;
  fish->step_ind = width * INIT_PERCENTAGE;
  fish->step_vol = fish->step_ind;
  fish->info.weight = W_SCALE / 2;
  srand(time(NULL));
  for (int i = 0; i < DIM_COUNT; i++) {
    fish->info.positions[i] = rand_real(func->params.init_min, func->params.init_max);
  }
  fish->info.value = func->f(fish->info.positions, DIM_COUNT);
  fish->info.value_improvement = 0;
}

void individual_move(fish_t *fish) {
  double next_pos[DIM_COUNT];
  compute_next_position(fish, next_pos);
  double next_val = fish->func.f(next_pos, DIM_COUNT);

  // Checks if the new position is better than the current one
  if (next_val < fish->info.value) {
    // The new position is better, so the fish moves
    for (int i = 0; i < DIM_COUNT; i++) {
      fish->info.displacements[i] = next_pos[i] - fish->info.positions[i];
      fish->info.positions[i] = next_pos[i];
    }
    fish->info.value_improvement = abs(fish->info.value - next_val);
    fish->info.value = next_val;
  }
}

// Updates weight of each fish based on its value improvement and the maximum
void feeding_operator( fish_t *fish, fish_info_t *fishes, int n) {
  fish->info.weight = compute_weight(fish, fishes, n);
  // MICHELE: as it is now, the weight of a fish can never decrease
  // Should we use the next_position before being set to 0 if value is not improved?
}

// Computes a weighted average of individual movements
void collective_instinctive_move(fish_t* fish, fish_info_t* fishes, int n) {
    // MICHELE: right now things here are computed by each fish, because we such information, but
    // we could parallelize things with an MPI_Allreduce-MPI_SUM or something similar
    // Is it worth it? I don't know
    
    double sum_displacements[DIM_COUNT] = {0};
    double sum_values = 0;

    // Compute the sum of all displacements and the sum of all values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < DIM_COUNT; j++) {
            sum_displacements[j] += fishes[i].displacements[j];
        }
        sum_values += fishes[i].value;
    }

    // Compute the collective instinctive move
    for (int j = 0; j < DIM_COUNT; j++) {
        double move = sum_displacements[j] / sum_values;
        fish->info.positions[j] += move;
    }
}

// Computes baricenter and move fishes towards/away from it
void collective_volitive_move(fish_t* fish, fish_info_t* fishes, int n) {
    // MICHELE: right now things here are computed by each fish, because we have such information, but
    // we could parallelize things with an MPI_Allreduce-MPI_SUM or something similar
    double baricenter[DIM_COUNT] = {0};
    double sum_weights = 0;

    // Compute the baricenter and the sum of all weights
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < DIM_COUNT; j++) {
            baricenter[j] += fishes[i].positions[j] * fishes[i].weight;
        }
        sum_weights += fishes[i].weight;
    }

    // TODO
    // Here we need to know wether the delta value is positive or negative
}

// Decreases step_ind and step_vol. To be used every cycle
void decrease_step(fish_t *fish, int cycle) {
  fish->step_ind = decrease_linearly(fish->step_ind, INIT_PERCENTAGE, FINAL_PERCENTAGE, CYCLES_LIMIT);
  fish->step_vol = decrease_linearly(fish->step_vol, INIT_PERCENTAGE, FINAL_PERCENTAGE, CYCLES_LIMIT);
}

/******************************************************************************/
