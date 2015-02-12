#include <mathfunctions.h>

void FindInterval(double a, double b, double *x, int N, int *imin, int *imax)
{
  int i;
  *imax = -1;
  *imin =  N;

  double min_dx = x[1] - x[0];
  for (i = 2; i < N; i++) {
    double dx = x[i] - x[i-1];
    if (dx < min_dx) min_dx = dx;
  }
  double tol = 1e-10 * min_dx;

  for (i = 0; i < N; i++) {
    if (x[i] <= (b+tol)) *imax = i+1;
  }
  for (i = N-1; i > -1; i--) {
    if (x[i] >= (a-tol)) *imin = i;
  }

  return;
}
