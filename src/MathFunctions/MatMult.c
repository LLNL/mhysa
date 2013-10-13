#include <mathfunctions.h>

/* A = XY */

int MatMult(int N,double **A,double **X,double **Y)
{
  int i,j,k;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i][j] = 0.0;
      for (k = 0; k < N; k++) A[i][j] += X[i][k]*Y[k][j];
    }
  }
  return(0);
}

