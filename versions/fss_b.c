#include "fss_b.h"
#include <stdlib.h>
#include <float.h>
#include <math.h>

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

void individual_move(fish_t* const fish, struct setup_info_t* const setup) {
  // Update food amount to current position (since it changed from collective
  // movements)
  double curr_val = compute_amount_of_food(setup->func.f(fish->positions, DIM_COUNT));

  double next_pos[DIM_COUNT];
  compute_next_position(fish, &setup->func, setup->step_ind, next_pos);

  // By making comparisons on the amount of food available in a position instead
  // of on the value of the functions being considered, we are allowed to always
  // look for the smallest possible value
  double next_val = compute_amount_of_food(setup->func.f(next_pos, DIM_COUNT));

  fish->food_improvement = next_val - curr_val;
  // Checks if the new position is better than the current one
  if (next_val >= curr_val) {
    // The new position is better, so the fish moves
    for (int i = 0; i < DIM_COUNT; i++) {
      fish->displacements[i] = next_pos[i] - fish->positions[i];
      fish->positions[i] = next_pos[i];
    }
  } else {
    // The new position is worse, so the fish stays in the current position
    for (int i = 0; i < DIM_COUNT; i++) {
      fish->displacements[i] = 0;
   }
  }
}

void feeding_operator(fish_t* const fish, double max_food_improvement) {
  double new_weight = compute_weight(fish, max_food_improvement);
  fish->weight_improvement = new_weight - fish->weight;
  fish->weight = new_weight;
}

void collective_instinctive_move(
  fish_t* const fish, double** displacements, double* food_improvements, int n,
  struct setup_info_t* const setup
) {
  double sum_displacements[DIM_COUNT] = {0};
  double total_value_improvement = 0;

  // Compute the sum of all displacements and the sum of all values
  for (int i = 0; i < n; i++) {
    double food_improvement = fmax(food_improvements[i], 0.0);
      for (int j = 0; j < DIM_COUNT; j++) {
        sum_displacements[j] += displacements[i][j] * food_improvement;
      }
      total_value_improvement += food_improvement;
  }

  // Compute the collective instinctive move
  if(total_value_improvement != 0.0) {
    for (int j = 0; j < DIM_COUNT; j++) {
      double move = sum_displacements[j] / total_value_improvement;
      fish->positions[j] += move;
      // Check that position is within bounds
      if (fish->positions[j] > setup->func.params.search_space_max) {
        fish->positions[j] = setup->func.params.search_space_max;
      } else if (fish->positions[j] < setup->func.params.search_space_min) {
        fish->positions[j] = setup->func.params.search_space_min;
      }
    }  
  }
}

void collective_volitive_move(
  fish_t* const fish, const fish_t* const fishes, int n,
  struct setup_info_t* const setup
) {
  double baricenter[DIM_COUNT] = {0};
  double total_weight_improvement = 0;
  double total_weight = 0;
//
  // Compute the baricenter and the sum of all weight_improvements
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < DIM_COUNT; j++) {
      baricenter[j] += fishes[i].positions[j] * fishes[i].weight;
    }
    total_weight += fishes[i].weight;
    total_weight_improvement += fishes[i].weight_improvement;
  }
//
  for (int j = 0; j < DIM_COUNT; j++) {
    baricenter[j] /= total_weight;
  }
//
  // If weights increased, we need to compact the group so that fishes move
  // towards the baricenter
  int inc = -1;
  if (total_weight_improvement < 0) {
    inc = +1;
  }
//
  // Compute the difference vector and its magnitude
  // NOTE: magnited leads to smother movement of fishes, but the end results
  // after some iterations is the same as when magnitude is not involved
  double diff[DIM_COUNT];
  double magnitude = 0;
  for (int j = 0; j < DIM_COUNT; j++) {
    diff[j] = fish->positions[j] - baricenter[j];
    magnitude += diff[j] * diff[j];
  }
  magnitude = sqrt(magnitude);
//
  // Normalize the difference vector
  if (magnitude > 0) {
    for (int j = 0; j < DIM_COUNT; j++) {
      diff[j] /= magnitude;
    }
  }
//
  for (int j = 0; j < DIM_COUNT; j++) {
    fish->positions[j] += inc * setup->step_vol * rand_real(0, 1) * diff[j];
    // Check that position is within bounds
    if (fish->positions[j] > setup->func.params.search_space_max) {
      fish->positions[j] = setup->func.params.search_space_max;
    } else if (fish->positions[j] < setup->func.params.search_space_min) {
      fish->positions[j] = setup->func.params.search_space_min;
    }
  }
}

void decrease_step(struct setup_info_t* setup) {
  setup->step_ind_perc -= setup->step_ind_perc_dec;
  setup->step_ind = setup->search_space_width * setup->step_ind_perc;
  setup->step_vol_perc -= setup->step_vol_perc_dec;
  setup->step_vol = setup->search_space_width * setup->step_vol_perc;
}

/******************************************************************************/
