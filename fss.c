#include "fss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/******************************************************************************/
/*** UTILITY FUNCTIONS *********************************************************/
/******************************************************************************/
/**
 * Returns a random value between a and b (exclusive).
 */
double rand_real(double min, double max) {
  double div = RAND_MAX / (max - min);
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

void compute_next_position(const fish_t* const fish, double next_pos[]) {
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

double compute_weight(const fish_t* const fish, const fish_info_t* const fishes, int n) {
  double max = -DBL_MAX;
  for (int i = 0; i < n; i++) {
    if (fishes[i].value_improvement > max) {
      max = fishes[i].value_improvement;
    }
  }
  return fish->info.weight + fish->info.value_improvement / max;
}

double decrease_linearly(double value, double initial_value, double final_value, int max_iterations) {
  // initial_value and final_value are percentages of the entire values space
  double dec_perc = (initial_value - final_value) / max_iterations;
  return value * (1 - dec_perc);
}

/******************************************************************************/
/*** FISH API & FSS MOVEMENT **************************************************/
/******************************************************************************/

void init(fish_t* const fish, struct func_t* const func) {
  double width = func->params.search_space_max - func->params.search_space_min;
  fish->step_ind = width * INIT_PERCENTAGE;
  fish->step_vol = fish->step_ind;
  fish->info.weight = W_SCALE / 2;
  for (int i = 0; i < DIM_COUNT; i++) {
    fish->info.positions[i] = rand_real(func->params.init_min, func->params.init_max);
  }
  fish->info.value = func->f(fish->info.positions, DIM_COUNT);
  fish->info.value_improvement = 0;
  fish->info.weight_improvement = 0;
  fish->func = *func;
}

void individual_move(fish_t* const fish) {
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
    fish->info.value_improvement = fish->info.value - next_val;
    fish->info.value = next_val;
  }
}

// Updates weight of each fish based on its value improvement and the maximum
void feeding_operator(fish_t* const fish, const fish_info_t* const fishes, int n) {
  fish->info.weight = compute_weight(fish, fishes, n);
}

// Computes a weighted average of individual movements
void collective_instinctive_move(fish_t* const fish, const fish_info_t* const fishes, int n) {
    // MICHELE: right now things here are computed by each fish, because we such
    // information, but we could parallelize things with an MPI_Allreduce-MPI_SUM
    // or something similar. Is it worth it? I don't know
    
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
void collective_volitive_move(fish_t* const fish, const fish_info_t* const fishes, int n, int i) {
    // MICHELE: right now things here are computed by each fish, because we have
    // such information, but we could parallelize things with an
    // MPI_Allreduce-MPI_SUM or something similar
    double baricenter[DIM_COUNT] = {0};
    double sum_weights = 0;

    // Compute the baricenter and the sum of all weight_improvements
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < DIM_COUNT; j++) {
            baricenter[j] += fishes[i].positions[j] * fishes[i].weight;
        }
        sum_weights += fishes[i].weight_improvement;
    }

    // if weigts increased, we need to compact the group having them go towards the baricenter
    int inc = -1;
    if (sum_weights < 0) {
      inc = +1;
    }

    // Compute the difference vector and its magnitude
    double diff[DIM_COUNT];
    double magnitude = 0;
    for (int j = 0; j < DIM_COUNT; j++) {
        diff[j] = baricenter[j] - fish->info.positions[j];
        magnitude += diff[j] * diff[j];
    }
    magnitude = sqrt(magnitude);

    // Normalize the difference vector
    if (magnitude > 0) {
        for (int j = 0; j < DIM_COUNT; j++) {
            diff[j] /= magnitude;
        }
    }

    for (int j = 0; j < DIM_COUNT; j++) {
      fish->info.positions[j] += inc * (fish->step_vol * rand_real(0, 1) * diff[j]);
    }
}

// Decreases step_ind and step_vol. To be used every cycle
void decrease_step(fish_t* const fish, int cycle) {
  fish->step_ind = decrease_linearly(fish->step_ind, INIT_PERCENTAGE, FINAL_PERCENTAGE, CYCLES_LIMIT);
  fish->step_vol = decrease_linearly(fish->step_vol, INIT_PERCENTAGE, FINAL_PERCENTAGE, CYCLES_LIMIT);
}

/******************************************************************************/
