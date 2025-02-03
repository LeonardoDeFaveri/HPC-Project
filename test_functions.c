#include "test_functions.h"
#include <math.h>
#include <stdio.h>

#define PI 3.1415926535897932384626433832
#define SQUARE(x) ((x) * (x))

double rosenbrock(double *x, int n) {
  double result = 0;
  for (int i = 0; i < n - 1; i++) {
    result += 100 * SQUARE(x[i + 1] - SQUARE(x[i])) + SQUARE(1 - x[i]);
  }
  return result;
}

double rastrigin(double *x, int n) {
  double result = 10 * n;
  for (int i = 0; i < n; i++) {
    result += SQUARE(x[i]) - 10 * cos(2 * PI * x[i]);
  }
  return result;
}

double griewank(double *x, int n) {
  double result = 0;
  for (int i = 0; i < n; i++) {
    result += SQUARE(x[i]);
  }
  result /= 4000;
  result += 1;
  double prod = 1;
  for (int i = 0; i < n; i++)
  {
    prod *= cos(x[i] / sqrt(i));
  }
  result -= prod;
  return result;
}

double ackley(double *x, int n) {
  double result = exp(1) + 20;
  double sum = 0;
  for (int i = 0; i < n; i++) {
    sum += cos(2 * PI * x[i]);
  }
  result -= exp(sum / n);
  return result;
}

double schwefel(double *x, int n) {
  double result = 0;
  for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j <= i; i++) {
      sum += x[j];
    }
    result += SQUARE(sum);
  }
  return result;
}

double empty(double *x, int n) {
  return 0.0;
}

const struct func_t get_function(const enum func_name name) {
  struct func_param_t params;
  double (*f)(double *, int);

  switch (name) {
  case ROSENBROCK:
    params = (struct func_param_t){-30.0, 30.0, 15.0, 30.0};
    f = *rosenbrock;
    break;
  case RASTRIGIN:
    params = (struct func_param_t){-5.12, 5.12, 2.56, 5.12};
    f = *rastrigin;
    break;
  case GRIEWANK:
    params = (struct func_param_t){-600.0, 600.0, 300.0, 600.0};
    f = *griewank;
    break;
  case ACKLEY:
    params = (struct func_param_t){-32.0, 32.0, 16.0, 32.0};
    f = *ackley;
    break;
  case SCHWEFEL:
    params = (struct func_param_t){-100.0, 100.0, 50.0, 100.0};
    f = *schwefel;
    break;
  default:
    f = *empty;
  }
  return (struct func_t) {name, params, f};
}