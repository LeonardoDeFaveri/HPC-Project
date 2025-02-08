#ifndef FSS_OLD_H
#define FSS_OLD_H

#include "../test_functions.h"

#define CYCLES_LIMIT 10000
#define W_SCALE 5000.0
/**
 * Test functions work on multidimensional inputs. This is the number of
 * dimentions.
 */
#define DIM_COUNT 30

/**
 * As the original authors did in the paper, `step_ind` and `step_vol` are set
 * as percentage of the actual search space. For initial values, percentages are
 * `10%` (`0.1`), `1%` (`0.01`) and `0.1%` (`0.001`).
 */
#define IND_INIT_PERCENTAGE 0.1
/**
 * Final values of `step_ind` and `step_vol` are set as percentages of the actual
 * search space. These percentages are `0.1%` (`0.001`), `0.01` (`0.0001`) and
 * `0.001` (`0.00001`).
 */
#define IND_FINAL_PERCENTAGE 0.001
#define VOL_INIT_PERCENTAGE 0.01
#define VOL_FINAL_PERCENTAGE 0.0001

struct fish_t {
  /**
   * Current position across all dimensions.
   */
  double positions[DIM_COUNT];
  /**
   * How has the position of the fish changed to get to the current one?
   */
  double displacements[DIM_COUNT];
  /**
   * How much weight does this fish have when deciding the collective displacement?
   */
  double weight;
  /**
   * How much weight changed after feeding ioperator.
   */
  double weight_improvement;
  /**
   * After a dispament, how much has the amount of food improved?
   */
  double food_improvement;
};
typedef struct fish_t fish_t;

struct setup_info_t {
  /**
   * Test function.
   */
  struct func_t func;
  /**
   * How much distance there is between the maximum allowed value in the search
   * space and the minimum one.
   */
  double search_space_width;
  /**
   * What percentage of the `search_space_width` is used at each step. 
   */
  double step_ind_perc;
  /**
   * How is `step_ind_perc` reduced at each step?
   */
  double step_ind_perc_dec;
  /**
   * How much displacement does the fish have when moving individually?
   */
  double step_ind;
  /**
   * What percentage of the `search_space_width` is used at each step. 
   */
  double step_vol_perc;
  /**
   * How is `step_ind_perc` reduced at each step?
   */
  double step_vol_perc_dec;
  /**
   * How much displacement does the fish have when moving together with the
   * other fishes?
   */
  double step_vol;
};

/**
 * Initializes an experiment setup using function `f`.
 */
void init_setup(struct setup_info_t* setup, const struct func_t* const f);
/**
 * Given an experiment setup initializes fish properties.
 * Fish initial position is randomly generated.
 */
void init(fish_t* const fish, const struct setup_info_t* const setup);

/******************************************************************************/
// FSS operations
void individual_move(fish_t* const fish, struct setup_info_t* const setup);
void feeding_operator(fish_t* const fish, double max_food_improvement);
void collective_instinctive_move(fish_t* const local_fishes, int local_count, struct setup_info_t* const setup);
void collective_volitive_move(
  fish_t* const local_fishes,   // array of local fishes
  int local_count,              // number of local fishes
  struct setup_info_t* const setup
);
void decrease_step(struct setup_info_t* setup);
/******************************************************************************/

#endif