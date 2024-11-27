#ifndef FSS_H
#define FSS_H

#define CYCLES_LIMIT 10000
#define W_SCALE 5000
/**
 * Test functions work on multidimensional inputs. This is the number of
 * dimentions.
 */
#define DIM_COUNT 5

struct fish_t {
  double step_ind;
  double step_vol;
  double func_val;
  double weight;
};
typedef struct fish_t fish_t;

#endif