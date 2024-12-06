#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

#include <stdint.h>

/**
 * Parameters required to set up a test environment for a function.
 */
struct func_param_t {
  double search_space_min;
  double search_space_max;
  double init_min;
  double init_max;
  double optima;
};

/**
 * Names of all supported functions.
 */
enum func_name {
  ROSENBROCK,
  RASTRIGIN,
  GRIEWANK,
  ACKLEY,
  SCHWEFEL,
  NO_FUNCTION,
};

/**
 * Bundles together a function with its parameters and its actual body.
 */
struct func_t {
  enum func_name name;
  struct func_param_t params;
  /**
   * Actual implementation of the function. Takes as input an array of doubles
   * and its size. This is to allow for multidimensional inputs.
   */
  double (*f)(double *, int);
};

double rosenbrock(double *x, int n);
double rastrigin(double *x, int n);
double griewank(double *x, int n);
double ackley(double *x, int n);
double schwefel(double *x, int n);
double empty(double *x, int n);

/**
 * Given a name, returns the associated function with the parameters
 * necessary for the test environment. If name is `NO_FUNCTION` or anything not
 * recognizable as the other variants, the result is a function that does nothing
 * if called and whose params have not been initialized, so they should not be
 * accessed.
 */
const struct func_t get_function(const enum func_name name);

#endif