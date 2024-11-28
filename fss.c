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
#define direction() (-1 * rand() % 2)

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

void collective_instinctive_move(fish_t *fish, fish_info_t *fishes, int n) {
  fish->info.weight = compute_weight(fish, fishes, n);
}

/******************************************************************************/
