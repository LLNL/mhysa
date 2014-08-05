/* Basic functions */
#include <math.h>

/* function to calculate the grid points corresponding to
 * a given interval */
void FindInterval      (double,double,double*,int,int*,int*);

/* minimum of two numbers */
#define min(a,b) ((a)<(b)?(a):(b))
/* maximum of two numbers */
#define max(a,b) ((a)>(b)?(a):(b))

/* minimum of three numbers */
#define min3(a,b,c) min(min((a),(b)),min((b),(c)))
/* maximum of three numbers */
#define max3(a,b,c) max(max((a),(b)),max((b),(c)))

/* absolute value */
#define absolute(a) ((a)<0?-(a):(a))

/* raise to a power: y = x^a */
#define raiseto(x,a) (exp((a)*log(x)))

/* return the sign of the argument */
#define sign(a) ((a)<0?-1.0:1.0)

/* Matrix-Matrix multiplication
 * N    : size of the matrices (NxN) (int)
 * A    : product matrix (double[])
 * X,Y  : input matrices (A = XY) (double[])
 *
 * Note: A, X, Y are 1D arrays of size N*N
 * saved in row-major format 
*/
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

/* Matrix-Matrix multiplication: loop unrolled for N=3 */
#define MatMult3(N,A,X,Y) \
  { \
    A[0] = X[0]*Y[0] + X[1]*Y[3] + X[2]*Y[6]; \
    A[1] = X[0]*Y[1] + X[1]*Y[4] + X[2]*Y[7]; \
    A[2] = X[0]*Y[2] + X[1]*Y[5] + X[2]*Y[8]; \
    A[3] = X[3]*Y[0] + X[4]*Y[3] + X[5]*Y[6]; \
    A[4] = X[3]*Y[1] + X[4]*Y[4] + X[5]*Y[7]; \
    A[5] = X[3]*Y[2] + X[4]*Y[5] + X[5]*Y[8]; \
    A[6] = X[6]*Y[0] + X[7]*Y[3] + X[8]*Y[6]; \
    A[7] = X[6]*Y[1] + X[7]*Y[4] + X[8]*Y[7]; \
    A[8] = X[6]*Y[2] + X[7]*Y[5] + X[8]*Y[8]; \
  }

/* Matrix-Matrix multiplication: loop unrolled for N=4 */
#define MatMult4(N,A,X,Y) \
  { \
    A[0]  =  X[0]*Y[0]  +  X[1]*Y[4]  +  X[2]*Y[8]  +  X[3]*Y[12]; \
    A[1]  =  X[0]*Y[1]  +  X[1]*Y[5]  +  X[2]*Y[9]  +  X[3]*Y[13]; \
    A[2]  =  X[0]*Y[2]  +  X[1]*Y[6]  +  X[2]*Y[10] +  X[3]*Y[14]; \
    A[3]  =  X[0]*Y[3]  +  X[1]*Y[7]  +  X[2]*Y[11] +  X[3]*Y[15]; \
    A[4]  =  X[4]*Y[0]  +  X[5]*Y[4]  +  X[6]*Y[8]  +  X[7]*Y[12]; \
    A[5]  =  X[4]*Y[1]  +  X[5]*Y[5]  +  X[6]*Y[9]  +  X[7]*Y[13]; \
    A[6]  =  X[4]*Y[2]  +  X[5]*Y[6]  +  X[6]*Y[10] +  X[7]*Y[14]; \
    A[7]  =  X[4]*Y[3]  +  X[5]*Y[7]  +  X[6]*Y[11] +  X[7]*Y[15]; \
    A[8]  =  X[8]*Y[0]  +  X[9]*Y[4]  +  X[10]*Y[8]  + X[11]*Y[12];\
    A[9]  =  X[8]*Y[1]  +  X[9]*Y[5]  +  X[10]*Y[9]  + X[11]*Y[13];\
    A[10] =  X[8]*Y[2]  +  X[9]*Y[6]  +  X[10]*Y[10] + X[11]*Y[14];\
    A[11] =  X[8]*Y[3]  +  X[9]*Y[7]  +  X[10]*Y[11] + X[11]*Y[15];\
    A[12] =  X[12]*Y[0] +  X[13]*Y[4] +  X[14]*Y[8]  + X[15]*Y[12];\
    A[13] =  X[12]*Y[1] +  X[13]*Y[5] +  X[14]*Y[9]  + X[15]*Y[13];\
    A[14] =  X[12]*Y[2] +  X[13]*Y[6] +  X[14]*Y[10] + X[15]*Y[14];\
    A[15] =  X[12]*Y[3] +  X[13]*Y[7] +  X[14]*Y[11] + X[15]*Y[15];\
  }

/* Matrix-Matrix multiplication: loop unrolled for N=5 */
#define MatMult5(N,A,X,Y) \
  { \
    A[0]  =  X[0]*Y[0]  +  X[1]*Y[5]  +  X[2]*Y[10]  +  X[3]*Y[15]  +  X[4]*Y[20]; \
    A[1]  =  X[0]*Y[1]  +  X[1]*Y[6]  +  X[2]*Y[11]  +  X[3]*Y[16]  +  X[4]*Y[21]; \
    A[2]  =  X[0]*Y[2]  +  X[1]*Y[7]  +  X[2]*Y[12]  +  X[3]*Y[17]  +  X[4]*Y[22]; \
    A[3]  =  X[0]*Y[3]  +  X[1]*Y[8]  +  X[2]*Y[13]  +  X[3]*Y[18]  +  X[4]*Y[23]; \
    A[4]  =  X[0]*Y[4]  +  X[1]*Y[9]  +  X[2]*Y[14]  +  X[3]*Y[19]  +  X[4]*Y[24]; \
    A[5]  =  X[5]*Y[0]  +  X[6]*Y[5]  +  X[7]*Y[10]  +  X[8]*Y[15]  +  X[9]*Y[20]; \
    A[6]  =  X[5]*Y[1]  +  X[6]*Y[6]  +  X[7]*Y[11]  +  X[8]*Y[16]  +  X[9]*Y[21]; \
    A[7]  =  X[5]*Y[2]  +  X[6]*Y[7]  +  X[7]*Y[12]  +  X[8]*Y[17]  +  X[9]*Y[22]; \
    A[8]  =  X[5]*Y[3]  +  X[6]*Y[8]  +  X[7]*Y[13]  +  X[8]*Y[18]  +  X[9]*Y[23]; \
    A[9]  =  X[5]*Y[4]  +  X[6]*Y[9]  +  X[7]*Y[14]  +  X[8]*Y[19]  +  X[9]*Y[24]; \
    A[10] =  X[10]*Y[0] +  X[11]*Y[5] +  X[12]*Y[10] +  X[13]*Y[15] +  X[14]*Y[20];\
    A[11] =  X[10]*Y[1] +  X[11]*Y[6] +  X[12]*Y[11] +  X[13]*Y[16] +  X[14]*Y[21];\
    A[12] =  X[10]*Y[2] +  X[11]*Y[7] +  X[12]*Y[12] +  X[13]*Y[17] +  X[14]*Y[22];\
    A[13] =  X[10]*Y[3] +  X[11]*Y[8] +  X[12]*Y[13] +  X[13]*Y[18] +  X[14]*Y[23];\
    A[14] =  X[10]*Y[4] +  X[11]*Y[9] +  X[12]*Y[14] +  X[13]*Y[19] +  X[14]*Y[24];\
    A[15] =  X[15]*Y[0] +  X[16]*Y[5] +  X[17]*Y[10] +  X[18]*Y[15] +  X[19]*Y[20];\
    A[16] =  X[15]*Y[1] +  X[16]*Y[6] +  X[17]*Y[11] +  X[18]*Y[16] +  X[19]*Y[21];\
    A[17] =  X[15]*Y[2] +  X[16]*Y[7] +  X[17]*Y[12] +  X[18]*Y[17] +  X[19]*Y[22];\
    A[18] =  X[15]*Y[3] +  X[16]*Y[8] +  X[17]*Y[13] +  X[18]*Y[18] +  X[19]*Y[23];\
    A[19] =  X[15]*Y[4] +  X[16]*Y[9] +  X[17]*Y[14] +  X[18]*Y[19] +  X[19]*Y[24];\
    A[20] =  X[20]*Y[0] +  X[21]*Y[5] +  X[22]*Y[10] +  X[23]*Y[15] +  X[24]*Y[20];\
    A[21] =  X[20]*Y[1] +  X[21]*Y[6] +  X[22]*Y[11] +  X[23]*Y[16] +  X[24]*Y[21];\
    A[22] =  X[20]*Y[2] +  X[21]*Y[7] +  X[22]*Y[12] +  X[23]*Y[17] +  X[24]*Y[22];\
    A[23] =  X[20]*Y[3] +  X[21]*Y[8] +  X[22]*Y[13] +  X[23]*Y[18] +  X[24]*Y[23];\
    A[24] =  X[20]*Y[4] +  X[21]*Y[9] +  X[22]*Y[14] +  X[23]*Y[19] +  X[24]*Y[24];\
  }

/* Matrix-Vector multiplication
 * N    : size of the matrix and vectors (int)
 * A    : product matrix (double[])
 * x    : input vector (double[])
 * y    : output vector (y = Ax) (double[])
 *
 * Note: A is 1D arrays of size N*N saved in row-major format,
 * x,y are 1D arrays of size N
*/
#define MatVecMult(N,y,A,x) \
  { \
    int i,j; \
    for (i = 0; i < (N); i++) { \
      y[i] = 0; \
      for (j = 0; j < (N); j++) y[i] += A[i*(N)+j]*x[j]; \
    } \
  }

/* Matrix-Vector multiplication: loop-unrolled for N=3 */
#define MatVecMult3(N,y,A,x) \
  { \
    y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];\
    y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];\
    y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];\
  }

/* Matrix-Vector multiplication: loop-unrolled for N=4 */
#define MatVecMult4(N,y,A,x) \
  { \
    y[0] = A[0]*x[0]  +  A[1]*x[1]  +  A[2]*x[2]  +  A[3]*x[3]; \
    y[1] = A[4]*x[0]  +  A[5]*x[1]  +  A[6]*x[2]  +  A[7]*x[3]; \
    y[2] = A[8]*x[0]  +  A[9]*x[1]  +  A[10]*x[2] +  A[11]*x[3];\
    y[3] = A[12]*x[0] +  A[13]*x[1] +  A[14]*x[2] +  A[15]*x[3];\
  }

/* Matrix-Vector multiplication: loop-unrolled for N=5 */
#define MatVecMult5(N,y,A,x) \
  { \
    y[0] = A[0]*x[0]  +  A[1]*x[1]  +  A[2]*x[2]  +  A[3]*x[3]  +  A[4]*x[4]; \
    y[1] = A[5]*x[0]  +  A[6]*x[1]  +  A[7]*x[2]  +  A[8]*x[3]  +  A[9]*x[4]; \
    y[2] = A[10]*x[0] + A[11]*x[1]  +  A[12]*x[2] +  A[13]*x[3] +  A[14]*x[4];\
    y[3] = A[15]*x[0] +  A[16]*x[1] +  A[17]*x[2] +  A[18]*x[3] +  A[19]*x[4];\
    y[4] = A[20]*x[0] +  A[21]*x[1] +  A[22]*x[2] +  A[23]*x[3] +  A[24]*x[4];\
  }
