#include "fss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

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
    // TODO place these constraints on every move?
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
    if (fishes[i].food_improvement > max) {
      max = fishes[i].food_improvement;
    }
  }

  // Comute the new weight, make sure its between bounds
  double new_weight = fish->info.weight;
  if (max != 0) {
    new_weight += fish->info.food_improvement / max;
    if (new_weight < 1)
      new_weight = 1;
    else if (new_weight > W_SCALE)
      new_weight = W_SCALE;
  }
  return new_weight;
}

double decrease_linearly(double value, double initial_value, double final_value, int max_iterations) {
  // initial_value and final_value are percentages of the entire values space
  double dec_perc = (initial_value - final_value) / max_iterations;
  return value * (1 - dec_perc);
}

/// The closer `value` is to `optimum_value`, the higher the returned value.
double compute_amount_of_food(double value, double optimum_value) {
  return DBL_MAX - abs(value - optimum_value);
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
  fish->info.food_amount = compute_amount_of_food(
    func->f(fish->info.positions, DIM_COUNT),
    fish->func.params.optimum
  );
  fish->info.food_improvement = 0;
  fish->info.weight_improvement = 0;
  fish->func = *func;
}

void individual_move(fish_t* const fish) {
  double next_pos[DIM_COUNT];
  compute_next_position(fish, next_pos);
  // By making comparisons on the amount of food available in a position instead
  // of on the value of the functions being considered, we are allowed to always
  // look for the highest possible value
  double next_val = compute_amount_of_food(
    fish->func.f(next_pos, DIM_COUNT),
    fish->func.params.optimum
  );

  // Checks if the new position is better than the current one
  if (next_val > fish->info.food_amount) {
    // The new position is better, so the fish moves
    for (int i = 0; i < DIM_COUNT; i++) {
      fish->info.displacements[i] = next_pos[i] - fish->info.positions[i];
      fish->info.positions[i] = next_pos[i];
    }
    fish->info.food_amount = next_val;
  } else {
    // The new position is worse, so the fish stays in the current position
    for (int i = 0; i < DIM_COUNT; i++) {
      fish->info.displacements[i] = 0;
    }
  }
  fish->info.food_improvement = next_val - fish->info.food_amount;
}

// Updates weight of each fish based on its value improvement and the maximum
void feeding_operator(fish_t* const fish, const fish_info_t* const fishes, int n) {
  double new_weight = compute_weight(fish, fishes, n);
  fish->info.weight_improvement = new_weight - fish->info.weight;
  fish->info.weight = new_weight;
}

// Computes a weighted average of individual movements
void collective_instinctive_move(fish_t* const fish, const fish_info_t* const fishes, int n) {
    // MICHELE: right now things here are computed by each fish, because we such
    // information, but we could parallelize things with an MPI_Allreduce-MPI_SUM
    // or something similar. Is it worth it? I don't know
    
    double sum_displacements[DIM_COUNT] = {0};
    double total_value_improvement = 0;

    // Compute the sum of all displacements and the sum of all values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < DIM_COUNT; j++) {
            sum_displacements[j] += fishes[i].displacements[j] * fishes[i].food_improvement;
        }
        total_value_improvement += fishes[i].food_improvement;
    }

    // Compute the collective instinctive move
    for (int j = 0; j < DIM_COUNT; j++) {
      if(total_value_improvement != 0.0) {
        double move = sum_displacements[j] / total_value_improvement;
        fish->info.positions[j] += move;
        // Check that position is within bounds
        if (fish->info.positions[j] > fish->func.params.search_space_max) {
          fish->info.positions[j] = fish->func.params.search_space_max;
        } else if (fish->info.positions[j] < fish->func.params.search_space_min) {
          fish->info.positions[j] = fish->func.params.search_space_min;
        }
      }  
    }
}

// Computes baricenter and move fishes towards/away from it
void collective_volitive_move(fish_t* const fish, const fish_info_t* const fishes, int n, int i) {
    // MICHELE: right now things here are computed by each fish, because we have
    // such information, but we could parallelize things with an
    // MPI_Allreduce-MPI_SUM or something similar
    double baricenter[DIM_COUNT] = {0};
    double total_positions[DIM_COUNT] = {0};
    double total_weight_improvement = 0;
    double total_weight = 0;

    // Compute the baricenter and the sum of all weight_improvements
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < DIM_COUNT; j++) {
        baricenter[j] += fishes[i].positions[j] * fishes[i].weight;
        total_positions[j] += fishes[i].positions[j];
      }
      total_weight += fishes[i].weight;
      total_weight_improvement += fishes[i].weight_improvement;
    }

    for (int j = 0; j < DIM_COUNT; j++) {
      baricenter[j] /= total_positions[j];
      // Normalize the baricenter by the total weight
      baricenter[j] /= total_weight;
    }

    // if weigts increased, we need to compact the group having them go towards the baricenter
    int inc = -1;
    if (total_weight_improvement < 0) {
      inc = +1;
    }

    // Compute the difference vector and its magnitude
    double diff[DIM_COUNT];
    double magnitude = 0;
    for (int j = 0; j < DIM_COUNT; j++) {
        //diff[j] = baricenter[j] - fish->info.positions[j];
        diff[j] = fish->info.positions[j] - baricenter[j];
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
      // Check that position is within bounds
      if (fish->info.positions[j] > fish->func.params.search_space_max) {
        fish->info.positions[j] = fish->func.params.search_space_max;
      } else if (fish->info.positions[j] < fish->func.params.search_space_min) {
        fish->info.positions[j] = fish->func.params.search_space_min;
      }
    }
}

// Decreases step_ind and step_vol. To be used every cycle
void decrease_step(fish_t* const fish, int cycle) {
  fish->step_ind = decrease_linearly(fish->step_ind, INIT_PERCENTAGE, FINAL_PERCENTAGE, CYCLES_LIMIT);
  fish->step_vol = decrease_linearly(fish->step_vol, INIT_PERCENTAGE, FINAL_PERCENTAGE, CYCLES_LIMIT);
}

/******************************************************************************/
