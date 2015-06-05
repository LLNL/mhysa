/*! @file matops.h
    @brief Contains macros and function definitions for common matrix operations.
    @author Debojyoti Ghosh
 */

#include <stdio.h>

/*! \def _MatrixZero_
 * Set all elements of the square matrix \a A of size \a N to 0, stored as 
 * a 1D array in row-major format.
*/
#define _MatrixZero_(A,N) \
  { \
    int arraycounter; \
    for (arraycounter=0; arraycounter<(N)*(N); arraycounter++) *((A)+arraycounter) = 0.0; \
  }

/*! \def _MatrixMultiply_
 * \a C = \a A \a B, where \a A, \a B, \a C are square matrices of size \a N,
 * stored as 1D arrays in row major format.
*/
#define _MatrixMultiply_(A,B,C,N) \
  { \
    int matopsi,matopsj,matopsk; \
    for (matopsi=0; matopsi<(N); matopsi++) \
      for (matopsj=0; matopsj<(N); matopsj++) { \
        *((C)+matopsi*(N)+matopsj) = 0; \
        for (matopsk=0; matopsk<(N); matopsk++) *((C)+matopsi*(N)+matopsj) += ((*((A)+matopsi*(N)+matopsk)) * (*((B)+matopsk*(N)+matopsj))); \
      } \
  }

/*! \def _MatrixMultiplyNonSquare_
 * \a C = \a A \a B where \a A is a \a NRowA X \a NColA matrix, 
 * \a B is a \a NColA X \a NColB matrix, and \a C is a \a NRowA X \a NColB
 * matrix. All matrices are stored as 1D arrays in row-major format.
*/
#define _MatrixMultiplyNonSquare_(A,B,C,NRowA,NColA,NColB) \
  { \
    int matopsi,matopsj,matopsk; \
    for (matopsi=0; matopsi<(NRowA); matopsi++) \
      for (matopsj=0; matopsj<(NColB); matopsj++) { \
        *((C)+matopsi*(NColB)+matopsj) = 0; \
        for (matopsk=0; matopsk<(NColA); matopsk++) *((C)+matopsi*(NColB)+matopsj) += ((*((A)+matopsi*(NColA)+matopsk)) * (*((B)+matopsk*(NColB)+matopsj))); \
      } \
  }

/*! \def _MatrixMultiplySubtract_
 * \a C = \a C - \a A \a B, where \a A, \a B, \a C are square matrices
 * of size \a N, saved as 1D arrays in row-major format.
*/
#define _MatrixMultiplySubtract_(C,A,B,N) \
  { \
    int matopsi,matopsj,matopsk; \
    for (matopsi=0; matopsi<(N); matopsi++) \
      for (matopsj=0; matopsj<(N); matopsj++) \
        for (matopsk=0; matopsk<(N); matopsk++) *((C)+matopsi*(N)+matopsj) -= ((*((A)+matopsi*(N)+matopsk)) * (*((B)+matopsk*(N)+matopsj))); \
  } \

/*! \def _MatVecMultiply_
 * \a y = \a A \a x, where \a x, \a y are vectors of size \a N, and
 * \a A is a square matrix of size \a N stored as a 1D array in row
 * major format.
*/
#define _MatVecMultiply_(A,x,y,N) \
  { \
    int matopsi,matopsj; \
    for (matopsi=0; matopsi<(N); matopsi++) { \
      *((y)+matopsi) = 0; \
      for (matopsj=0; matopsj<(N); matopsj++) *((y)+matopsi) += (*((A)+matopsi*N+matopsj) * *((x)+matopsj)); \
    } \
  }

/*! \def _MatVecMultiplySubtract_
 * \a y = \a y - \a A \a x, where \a x, \a y are vectors of size \a N, and
 * \a A is a square matrix of size \a N stored as a 1D array in row
 * major format.
*/
#define _MatVecMultiplySubtract_(y,A,x,N) \
  { \
    int matopsi,matopsj; \
    for (matopsi=0; matopsi<N; matopsi++) \
      for (matopsj=0; matopsj<(N); matopsj++) *((y)+matopsi) -= (*((A)+matopsi*N+matopsj) * *((x)+matopsj)); \
  }

/*! \def _MatrixInvert_
 * B =A^{-1} (Inverse of \a A), where \a A, \a B are square matrices of
 * size \a N stored as 1D arrays in row-major format. Uses LU decomposition.
 * If \a A is singular, \a B will contain NaN or Inf.
*/
#define _MatrixInvert_(A,B,N)  \
  { \
    int matopsi,matopsj,matopsk; \
    double matopsfactor, matopssum, matopsAc[(N)*(N)]; \
    \
    /* make a copy of A */ \
    for (matopsi=0; matopsi<(N)*(N); matopsi++) *(matopsAc+matopsi) = *((A)+matopsi); \
    \
    /* set B as the identity matrix */  \
    for (matopsi=0; matopsi<(N)*(N); matopsi++) *((B)+matopsi)             = 0.0; \
    for (matopsi=0; matopsi<(N)    ; matopsi++) *((B)+matopsi*(N)+matopsi) = 1.0; \
    \
    /* LU Decomposition - Forward Sweep */ \
    for (matopsi=0; matopsi<(N)-1; matopsi++) { \
      if (*(matopsAc+matopsi*(N)+matopsi) == 0) fprintf(stderr,"Error in MatrixInvert(): Matrix is singular.\n"); \
      for (matopsj=matopsi+1; matopsj<(N); matopsj++) { \
        matopsfactor = *(matopsAc+matopsj*(N)+matopsi) / *(matopsAc+matopsi*(N)+matopsi); \
        for (matopsk=matopsi+1; matopsk<(N); matopsk++) *(matopsAc+matopsj*(N)+matopsk) -= (matopsfactor * *(matopsAc +matopsi*(N)+matopsk)); \
        for (matopsk=0; matopsk<matopsj; matopsk++) *((B)+matopsj*(N)+matopsk) -= (matopsfactor * *((B)+matopsi*(N)+matopsk)); \
      } \
    } \
    \
    /* LU Decomposition - Backward Sweep */ \
    for (matopsi=(N)-1; matopsi>=0; matopsi--) { \
      for (matopsk=0; matopsk<(N); matopsk++) { \
        matopssum = 0.0; \
        for (matopsj=matopsi+1; matopsj<(N); matopsj++) matopssum += (*(matopsAc+matopsi*(N)+matopsj) * *((B)+matopsj*(N)+matopsk)); \
        *((B)+matopsi*(N)+matopsk) = (*((B)+matopsi*(N)+matopsk) - matopssum) / *(matopsAc+matopsi*(N)+matopsi); \
      } \
    } \
    \
    /* Done - B contains A^{-1} now */ \
  }
