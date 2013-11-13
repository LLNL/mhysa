/* Basic functions */
#include <math.h>

void FindInterval      (double,double,double*,int,int*,int*);

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

#define min3(a,b,c) min(min((a),(b)),min((b),(c)))
#define max3(a,b,c) max(max((a),(b)),max((b),(c)))

#define absolute(a) ((a)<0?-(a):(a))

#define raiseto(x,a) (exp((a)*log(x)))

#define sign(a) ((a)<0?-1.0:1.0)

#define MatMult(N,A,X,Y) \
  { \
    int i,j,k;  \
    for (i = 0; i < (N); i++) { \
      for (j = 0; j < (N); j++) { \
        A[i*(N)+j] = 0.0; \
        for (k = 0; k < (N); k++) A[i*(N)+j] += X[i*(N)+k]*Y[k*(N)+j]; \
      } \
    } \
  }

#define MatVecMult(N,y,A,x) \
  { \
    int i,j; \
    for (i = 0; i < (N); i++) { \
      y[i] = 0; \
      for (j = 0; j < (N); j++) y[i] += A[i*(N)+j]*x[j]; \
    } \
  }
