#include <mathfunctions.h>

/* A = XY */

int MatMult(int N,double *A,double *X,double *Y)
{
  int i,j,k;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i*N+j] = 0.0;
      for (k = 0; k < N; k++) A[i*N+j] += X[i*N+k]*Y[k*N+j];
    }
  }
  return(0);
}

