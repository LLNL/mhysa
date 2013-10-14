#include <mathfunctions.h>

/* y = Ax */

inline int MatVecMult(int N,double *y,double **A,double *x)
{
  int i,j;
  for (i = 0; i < N; i++) {
    y[i] = 0;
    for (j = 0; j < N; j++) y[i] += A[i][j]*x[j];
  }
  return(0);
}
