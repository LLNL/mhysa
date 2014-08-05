#include <stdio.h>
#include <stdlib.h>
#include <basic.h>

/* Returns the 1D index in the 1D array, given an n-D index in an n-D array
 * Input arguments:-
 *  N     : number of dimensions (int)
 *  imax  : array of size N, with the size of the n-D array in each dimension (int[])
 *          (not including ghost points)
 *  i     : the n-D index (integer array of size N) (int[])
 *  ghosts: number of ghost points (int)
 * Output:-
 *  index : the 1-D index (int)
*/
#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

/* same as _ArrayIndex1D_, loop-unrolled for 2D */
#define _ArrayIndex1D2_(N,imax,i,ghost,index) \
  { \
    index = i[1]+(ghost); \
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost))); \
  }

/* same as _ArrayIndex1D_, loop-unrolled for 3D */
#define _ArrayIndex1D3_(N,imax,i,ghost,index)  \
  { \
    index = i[2]+(ghost); \
    index = ((index*(imax[1]+2*(ghost))) + (i[1]+(ghost))); \
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost))); \
  }

/* Returns the 1D index in the 1D array, given an n-D index with an offset in an n-D array
 * Input arguments:-
 *  N     : number of dimensions (int)
 *  imax  : array of size N, with the size of the n-D array in each dimension (int[])
 *          (not including ghost points)
 *  i     : the n-D index (integer array of size N) (int[])
 *  offset: array of size N containing the offsets in each dimension (int[])
 *  ghosts: number of ghost points (int)
 * Output:-
 *  index : the 1-D index (int)
*/
#define _ArrayIndex1DWO_(N,imax,i,offset,ghost,index) \
  { \
    index = i[N-1]+(ghost)+ offset[N-1];\
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost)+offset[arraycounter]));\
    } \
  }

/* same as _ArrayIndex1DWO_, loop-unrolled for 2D */
#define _ArrayIndex1DWO2_(N,imax,i,offset,ghost,index) \
  { \
    index = i[1]+(ghost)+ offset[1];\
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost)+offset[0]));\
  }

/* same as _ArrayIndex1DWO_, loop-unrolled for 3D */
#define _ArrayIndex1DWO3_(N,imax,i,offset,ghost,index) \
  { \
    index = i[2]+(ghost)+ offset[2];\
    index = ((index*(imax[1]+2*(ghost))) + (i[1]+(ghost)+offset[1]));\
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost)+offset[0]));\
  }

/* Increment the n-D index by one 
 * Input:-
 *  N     : number of dimensions (int)
 *  imax  : array of size N, with the size of the n-D array in each dimension (int[])
 *          (not including ghost points)
 *  i     : the n-D index (gets incremented) (int[])
 * Output:-
 *  done  : set to 1 if i has reached its maximum value, else 0 (int)
*/
#define _ArrayIncrementIndex_(N,imax,i,done) \
  { \
    int arraycounter = 0; \
    while (arraycounter < (N)) { \
      if (i[arraycounter] == imax[arraycounter]-1) { \
        i[arraycounter] = 0; \
        arraycounter++; \
      } else { \
        i[arraycounter]++; \
        break; \
      } \
    } \
    if (arraycounter == (N)) done = 1; \
    else          done = 0; \
  }

/* Set all elements of an array (any datatype) to a constant value
 * Arguments:-
 *  x     : the array (float/double/int [])
 *  size  : the size of the array (int)
 *  value : the constant value (float/double/int)
*/
#define _ArraySetValue_(x,size,value)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);                       \
  }

/* Multiple all elements of an array by a constant value
 * Arguments:-
 *  x     : the array (float/double/int [])
 *  a     : the constant factor to multiply with (float/double/int)
 *  size  : size of the array (int)
*/
#define _ArrayScale1D_(x,a,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] *= a;                                   \
  }

/* Element-wise subtraction x = a - b
 * Arguments:
 *  a,b,x   : the arrays (x=a-b) (int/float/double [])
 *  size    : size of the arrays (int)
*/
#define _ArraySubtract1D_(x,a,b,size)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] - b[arraycounter];    \
  }

/* Element-wise addition x = a + b
 * Arguments:
 *  a,b,x   : the arrays (x=a+b) (int/float/double [])
 *  size    : size of the arrays (int)
*/
#define _ArrayAdd1D_(x,a,b,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] + b[arraycounter];    \
  }

/* Element-wise AXPY y = y + a*x
 * Arguments:
 *  x,y     : the arrays (y=y+a*x) (int/float/double [])
 *  a       : the constant value (int/float/double)
 *  size    : size of the arrays (int)
*/
#define _ArrayAXPY_(x,a,y,size)                                                                                     \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) y[arraycounter] += a*x[arraycounter];                   \
  }

/* Element-wise copy y = x
 * x,y    : the arrays (int/float/double/char [])
 * size   : size of the arrays (int)
*/
#define _ArrayCopy1D_(x,y,size) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = x[arraycounter]; \
  }

/* Same as _ArrayCopy1D_, loop-unrolled for size = 2 */
#define _ArrayCopy1D2_(x,y,size) \
  { \
    y[0] = x[0]; \
    y[1] = x[1]; \
  }

/* Same as _ArrayCopy1D_, loop-unrolled for size = 3 */
#define _ArrayCopy1D3_(x,y,size) \
  { \
    y[0] = x[0]; \
    y[1] = x[1]; \
    y[2] = x[2]; \
  }

#if !defined(INLINE)
# define INLINE inline
#endif

INLINE int    ArrayCopynD            (int,double*,double*,int*,int,int,int*,int);
INLINE double ArrayMaxnD             (int,int,int*,int,int*,double*);
INLINE double ArraySumSquarenD       (int,int,int*,int,int*,double*);
INLINE double ArraySumAbsnD          (int,int,int*,int,int*,double*);

/* Copy one n-D array to another n-D array
 * Arguments:-
 *  ndims       : number of dimensions (int)
 *  x           : copy-from array (double[])
 *  y           : copy-to   array (double[])
 *  dim         : integer array of size in each dimension (int[])
 *  g1          : number of ghost points in copy-from array x (int)
 *  g2          : number of ghost points in copy-to   array y (int)
 *  index       : pre-allocated integer array of size ndims (int[])
 *  nvars       : number of elements at one array location (int)
 *                (can be > 1 for systems of equations)
*/
INLINE int ArrayCopynD(int ndims,double *x,double *y,int *dim,int g1,int g2,int *index,int nvars)
{
  if (!y) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"y\" not allocated.\n");
    return(1);
  }
  if (!x) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"x\" not allocated.\n");
    return(1);
  }
  int done = 0;
  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p1, p2, n;
    _ArrayIndex1D_(ndims,dim,index,g1,p1); 
    _ArrayIndex1D_(ndims,dim,index,g2,p2);
    for (n=0; n<nvars; n++) y[p2*nvars+n] = x[p1*nvars+n];
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(0);
}

/* Returns the maximum magnitude element in an n-D array (useful for L_inf norm)
 * Arguments:-
 *  nvars       : number of elements at one array location (int)
 *                (can be > 1 for systems of equations)
 *  ndims       : number of dimensions (int)
 *  dim         : integer array of size in each dimension (int[])
 *  ghosts      : number of ghost points in the array x (int)
 *  index       : pre-allocated integer array of size ndims (int[])
 *  x           : the array (double[])
*/
INLINE double ArrayMaxnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; 
    for (v=0; v<nvars; v++) {
      double term = ( x[p*nvars+v]>0 ? x[p*nvars+v] : -x[p*nvars+v] );
      if (term > sum) sum = term;
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

/* Returns the sum-of-magnitudes of the elements in an n-D array (useful for L_1 norm)
 * Arguments:-
 *  nvars       : number of elements at one array location (int)
 *                (can be > 1 for systems of equations)
 *  ndims       : number of dimensions (int)
 *  dim         : integer array of size in each dimension (int[])
 *  ghosts      : number of ghost points in the array x (int)
 *  index       : pre-allocated integer array of size ndims (int[])
 *  x           : the array (double[])
*/
INLINE double ArraySumAbsnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += ( x[p*nvars+v]>0 ? x[p*nvars+v] : -x[p*nvars+v] );
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

/* Returns the sum-of-squares of the elements in an n-D array (useful for L_2 norm)
 * Arguments:-
 *  nvars       : number of elements at one array location (int)
 *                (can be > 1 for systems of equations)
 *  ndims       : number of dimensions (int)
 *  dim         : integer array of size in each dimension (int[])
 *  ghosts      : number of ghost points in the array x (int)
 *  index       : pre-allocated integer array of size ndims (int[])
 *  x           : the array (double[])
*/
INLINE double ArraySumSquarenD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += (x[p*nvars+v]*x[p*nvars+v]);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}
