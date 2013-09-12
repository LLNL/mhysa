#include <mathfunctions.h>

void FindInterval(double a, double b, double *x, int N, int *imin, int *imax)
{
  int i;
  *imax = -1;
  *imin =  N;

  for (i = 0; i < N; i++) {
    if (x[i] <= b) *imax = i+1;
  }
  for (i = N-1; i > -1; i--) {
    if (x[i] >= a) *imin = i;
  }

  return;
}
