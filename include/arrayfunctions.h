/*! @file arrayfunctions.h
    @brief Contains macros and function definitions for common array operations.
    @author Debojyoti Ghosh
 */

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>

/*! \def _ArrayIndexnD_ 
 * Returns the \a N -dimensional index \a i[\a N] in an \a N -dimensional array of size 
 * \a imax[\a N] with \a ghost number of ghost points, given the 
 * 1-dimensional \a index in the corresponding 1-dimensional array stored in
 * memory.
 * \sa #_ArrayIndex1D_
*/
#define _ArrayIndexnD_(N,index,imax,i,ghost)  \
  { \
    int arraycounter, term=1, index_copy=index; \
    for (arraycounter=0; arraycounter<N; arraycounter++) term *= (imax[arraycounter]+2*(ghost));\
    for (arraycounter=N-1; arraycounter>=0; arraycounter--) {\
      term /= imax[arraycounter]; \
      i[arraycounter] = index_copy/term; \
      index_copy -= i[arraycounter]*term; \
    } \
    for (arraycounter=0; arraycounter<N; arraycounter++) i[arraycounter] -= (ghost);\
  }

/*! \def _ArrayIndex1D_ 
 * Returns 1-dimensional \a index in the corresponding 1-dimensional array 
 * stored in memory, given the \a N -dimensional index \a i[\a N] in an 
 * \a N -dimensional array of size \a imax[\a N] with \a ghost number of 
 * ghost points.
 * \sa #_ArrayIndexnD_, #_ArrayIndex1D2_, #_ArrayIndex1D3_
*/
#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

/*! \def _ArrayIndex1D2_ 
 * Returns 1-dimensional \a index in the corresponding 1-dimensional array 
 * stored in memory, given the 2-dimensional index \a i[2] in an 
 * 2-dimensional array of size \a imax[2] with \a ghost number of 
 * ghost points. \a N is not used, but retained for uniformity.
 * \sa #_ArrayIndex1D_, #_ArrayIndex1D3_
*/
#define _ArrayIndex1D2_(N,imax,i,ghost,index) \
  { \
    index = i[1]+(ghost); \
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost))); \
  }

/*! \def _ArrayIndex1D3_ 
 * Returns 1-dimensional \a index in the corresponding 1-dimensional array 
 * stored in memory, given the 3-dimensional index \a i[3] in an 
 * 3-dimensional array of size \a imax[3] with \a ghost number of 
 * ghost points. \a N is not used, but retained for uniformity.
 * \sa #_ArrayIndex1D_, #_ArrayIndex1D2_
*/
#define _ArrayIndex1D3_(N,imax,i,ghost,index)  \
  { \
    index = i[2]+(ghost); \
    index = ((index*(imax[1]+2*(ghost))) + (i[1]+(ghost))); \
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost))); \
  }

/*! \def _ArrayIndex1DWO_ 
 * Returns 1-dimensional \a index in the corresponding 1-dimensional array 
 * stored in memory, given the \a N -dimensional index \a i[\a N] + an 
 * \a offset[\a N] in an \a N -dimensional array of size \a imax[\a N] 
 * with \a ghost number of ghost points.
 * \sa #_ArrayIndexnD_, #_ArrayIndex1D_, #_ArrayIndex1DWO2_, #_ArrayIndex1DWO3_
*/
#define _ArrayIndex1DWO_(N,imax,i,offset,ghost,index) \
  { \
    index = i[N-1]+(ghost)+ offset[N-1];\
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost)+offset[arraycounter]));\
    } \
  }

/*! \def _ArrayIndex1DWO2_ 
 * Returns 1-dimensional \a index in the corresponding 1-dimensional array 
 * stored in memory, given the 2-dimensional index \a i[2] + an 
 * \a offset[2] in an 2-dimensional array of size \a imax[2] 
 * with \a ghost number of ghost points. \a N is not used, but retained
 * for uniformity.
 * \sa #_ArrayIndex1DWO_, #_ArrayIndex1DWO3_
*/
#define _ArrayIndex1DWO2_(N,imax,i,offset,ghost,index) \
  { \
    index = i[1]+(ghost)+ offset[1];\
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost)+offset[0]));\
  }

/*! \def _ArrayIndex1DWO3_ 
 * Returns 1-dimensional \a index in the corresponding 1-dimensional array 
 * stored in memory, given the 3-dimensional index \a i[3] + an 
 * \a offset[3] in an 3-dimensional array of size \a imax[3] 
 * with \a ghost number of ghost points. \a N is not used, but retained
 * for uniformity.
 * \sa #_ArrayIndex1DWO_, #_ArrayIndex1DWO2_
*/
#define _ArrayIndex1DWO3_(N,imax,i,offset,ghost,index) \
  { \
    index = i[2]+(ghost)+ offset[2];\
    index = ((index*(imax[1]+2*(ghost))) + (i[1]+(ghost)+offset[1]));\
    index = ((index*(imax[0]+2*(ghost))) + (i[0]+(ghost)+offset[0]));\
  }

/*! \def _ArrayIncrementIndex_
 * Increments an \a N -dimensional index \a i[\a N] by one. If it reaches
 * the provided bounds \a imax[\a N], i.e., if \a i[c] = \a imax[c]-1 for all 
 * c = 0,...,\a N-1), then \a done = 1; else \a done = 0.
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

/*! \def _ArraySetValue_
 * Set all elements of a 1-dimensional array \a x (any datatype)
 * of length \a size to a scalar \a value
*/
#define _ArraySetValue_(x,size,value)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);                       \
  }

/*! \def _ArrayScale1D_ 
 * Multiply all elements of a 1-dimensional array \a x of length 
 * \a size by a scalar \a value
*/
#define _ArrayScale1D_(x,a,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] *= a;                                   \
  }

/*! \def _ArraySubtract1D_
 * Element-wise subtraction \a x = \a a - \a b, where \a x, \a a, 
 * and \a b are 1-dimensional arrays of length \a size
 * \sa #_ArrayAdd1D_, #_ArrayMultiply1D_
*/
#define _ArraySubtract1D_(x,a,b,size)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] - b[arraycounter];    \
  }

/*! \def _ArrayAdd1D_
 * Element-wise addition \a x = \a a + \a b, where \a x, \a a, 
 * and \a b are 1-dimensional arrays of length \a size
 * \sa #_ArraySubtract1D_, #_ArrayMultiply1D_
*/
#define _ArrayAdd1D_(x,a,b,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] + b[arraycounter];    \
  }

/*! \def _ArrayMultiply1D_
 * Element-wise multiplication \a x = \a a . \a b, where \a x, \a a, 
 * and \a b are 1-dimensional arrays of length \a size (a.k.a. dot product).
 * \sa #_ArraySubtract1D_, #_ArrayAdd1D_
*/
#define _ArrayMultiply1D_(x,a,b,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] * b[arraycounter];    \
  }

/*! \def _ArrayMultiply3Add1D_
 * Element-wise \a x =\a a . \a b + \a c . \a d + \a e . \a f, where \a x, \a a, 
 * \a b, \a c, \a d, \a e, and \a f, are 1-dimensional arrays of length \a size.
 * \sa #_ArrayMultiply1D_
*/
#define _ArrayMultiply3Add1D_(x,a,b,c,d,e,f,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) \
      x[arraycounter] = a[arraycounter]*b[arraycounter]+c[arraycounter]*d[arraycounter]+e[arraycounter]*f[arraycounter]; \
  }

/*! \def _ArrayConvexCombination1D_
 * Element-wise convex combination \a z = \a a . \a x + (1-\a a) . \a y, 
 * where \a a, \a x, \a y, \a z are 1-dimensional arrays of length \a size.
*/
#define _ArrayConvexCombination1D_(z,a,x,y,size)                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) \
      z[arraycounter] = a[arraycounter]*x[arraycounter]+(1.0-a[arraycounter])*y[arraycounter]; \
  }

/*! \def _ArrayAYPX_
 * Element-wise AYPX \a y = \a a \a y + \a x, 
 * where \a a is a scalar, and \a x, \a y, \a z are 
 * 1-dimensional arrays of length \a size.
 * \sa #_ArrayAXPY_, #_ArrayAXBY_, #_ArrayAXBYCZ_, #_ArrayScaledAXPY_
*/
#define _ArrayAYPX_(x,a,y,size)                                                                                     \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) \
      y[arraycounter] = a*y[arraycounter] + x[arraycounter];\
  }

/*! \def _ArrayAXPY_
 * Element-wise AXPY \a y = \a a \a x + \a y, 
 * where \a a is a scalar, and \a x, \a y, \a z are 
 * 1-dimensional arrays of length \a size.
 * \sa #_ArrayAYPX_, #_ArrayAXBY_, #_ArrayAXBYCZ_, #_ArrayScaledAXPY_
*/
#define _ArrayAXPY_(x,a,y,size)                                                                                     \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) y[arraycounter] += a*x[arraycounter];                   \
  }

/*! \def _ArrayAXBY_
 * Element-wise AXPY \a z = \a a \a x + \a b \a y, 
 * where \a a , \a b are scalars, and \a x, \a y, \a z are 
 * 1-dimensional arrays of length \a size.
 * \sa #_ArrayAYPX_, #_ArrayAXPY_, #_ArrayAXBYCZ_, #_ArrayScaledAXPY_
*/
#define _ArrayAXBY_(z,a,x,b,y,size)                                                                                 \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) z[arraycounter] = a*x[arraycounter]+b*y[arraycounter];  \
  }

/*! \def _ArrayAXBYCZ_
 * Element-wise AXPY \a w = \a a \a x + \a b \a y + \a c \a z, 
 * where \a a, \a b, \a c are scalars, and \a x, \a y, \a z, \a w
 * are 1-dimensional arrays of length \a size.
 * \sa #_ArrayAYPX_, #_ArrayAXPY_, #_ArrayScaledAXPY_
*/
#define _ArrayAXBYCZ_(w,a,x,b,y,c,z,size)                                                                                 \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) w[arraycounter] = a*x[arraycounter]+b*y[arraycounter]+c*z[arraycounter];  \
  }

/*! \def _ArrayScaledAXPY_
 *Element-wise Scaled AXPY \a y = \a e * (\a y + \a a * \a x), 
 * where \a a, \a e are scalars, and \a x, \a y
 * are 1-dimensional arrays of length \a size.
 * \sa #_ArrayAYPX_, #_ArrayAXPY_, #_ArrayAXBYCZ_
*/
#define _ArrayScaledAXPY_(x,a,e,y,size)                                                                                     \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) \
      y[arraycounter] = e*(y[arraycounter]+a*x[arraycounter]);                   \
  }

/*! \def _ArrayCopy1D_
 * Element-wise copy \a y = \a x, 
 * where \a x, \a y are 1-dimensional arrays of length \a size.
 * \sa #_ArrayCopy1D2_, #_ArrayCopy1D3_, #_ArrayScaleCopy1D_, #_ArrayAddCopy1D_
*/
#define _ArrayCopy1D_(x,y,size) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = x[arraycounter]; \
  }

/*! \def _ArrayCopy1D2_
 * Element-wise copy \a y = \a x, 
 * where \a x, \a y are 1-dimensional arrays of length 2. 
 * \a size is not used but retained for uniformity.
 * \sa #_ArrayCopy1D_, #_ArrayCopy1D3_
*/
#define _ArrayCopy1D2_(x,y,size) \
  { \
    y[0] = x[0]; \
    y[1] = x[1]; \
  }

/*! \def _ArrayCopy1D3_
 * Element-wise copy \a y = \a x, 
 * where \a x, \a y are 1-dimensional arrays of length 3. 
 * \a size is not used but retained for uniformity.
 * \sa #_ArrayCopy1D_, #_ArrayCopy1D2_
*/
#define _ArrayCopy1D3_(x,y,size) \
  { \
    y[0] = x[0]; \
    y[1] = x[1]; \
    y[2] = x[2]; \
  }

/*! \def _ArrayScaleCopy1D_
 * Element-wise scale and copy \a y = \a a * \a x, where
 * \a is a scalar, and \a x, \a y are one-dimensional arrays
 * of length \a size.
 * \sa #_ArrayCopy1D_, #_ArrayAddCopy1D_
*/
#define _ArrayScaleCopy1D_(x,a,y,size) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = a*x[arraycounter]; \
  }

/*! \def _ArrayAddCopy1D_ 
 * Element-wise add and copy \a y = \a a + \a x, 
 * where \a x,\a y are one-dimensional arrays of length \a size,
 * and \a a is a scalar.
 * \sa #_ArrayCopy1D_, #_ArrayScaleCopy1D_
*/
#define _ArrayAddCopy1D_(x,a,y,size) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = a+x[arraycounter]; \
  }

/*! \def _ArrayProduct1D_
 * \a p is the product of all the elements of the one-dimensional 
 * array \a x of length \a size. */
#define _ArrayProduct1D_(x,size,p) \
  { \
    int arraycounter = 0; p = 1; \
    for (arraycounter=0; arraycounter<size; arraycounter++) p *= x[arraycounter]; \
  }

/*! \def _ArrayBlockMultiply_
  Given two arrays: \a x of size \a n*bs, and \a a of size \a n, this macro
  implements: \a x[i][j] *= \a a[i] where \a i = 1,...,\a n, j = 1,...,\a bs,
  and \a x is stored as a 1D array in row-major format, i.e., \a x[i][j] = \a x[i*bs+j].
*/
#define _ArrayBlockMultiply_(x,a,n,bs) \
  { \
    int arraycounter1,arraycounter2; \
    for (arraycounter1=0; arraycounter1<n; arraycounter1++) { \
      for (arraycounter2=0; arraycounter2<bs; arraycounter2++) { \
        x[bs*arraycounter1+arraycounter2] *= a[arraycounter1]; \
      } \
    }\
  }

#if !defined(INLINE)
# define INLINE inline
#endif

INLINE int    ArrayCopynD            (int,double*,double*,int*,int,int,int*,int);
INLINE double ArrayMaxnD             (int,int,int*,int,int*,double*);
INLINE double ArraySumSquarenD       (int,int,int*,int,int*,double*);
INLINE double ArraySumAbsnD          (int,int,int*,int,int*,double*);

/*! Copy one n-D array to another n-D array (both of which are stored in memory as 1D arrays) */
INLINE int ArrayCopynD(int    ndims,  /*!< number of dimensions */
                       double *x,     /*!< copy-from array */
                       double *y,     /*!< copy-to   array */
                       int    *dim,   /*!< integer array of size in each dimension */
                       int    g1,     /*!< number of ghost points in copy-from array x */
                       int    g2,     /*!< number of ghost points in copy-to   array y */
                       int    *index, /*!< pre-allocated (by the calling function) integer array of size ndims */
                       int    nvars   /*!< number of elements at one array location,
                                           can be > 1 for systems of equations)*/
                      )
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
    int p1, p2;
    _ArrayIndex1D_(ndims,dim,index,g1,p1); 
    _ArrayIndex1D_(ndims,dim,index,g2,p2);
    _ArrayCopy1D_((x+p1*nvars),(y+p2*nvars),nvars);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(0);
}

/*! Returns the maximum magnitude element in an n-D array (useful for L_inf norm) */
INLINE double ArrayMaxnD(int    nvars,  /*!< number of elements at one array location,
                                             can be > 1 for systems of equations */
                         int    ndims,  /*!< number of dimensions */
                         int    *dim,   /*!< integer array containing the size of x in each dimension*/
                         int    ghosts, /*!< number of ghost points in the array x*/
                         int    *index, /*!< pre-allocated (by the calling function) integer array of size ndims */
                         double *x      /*!< the array*/
                        )
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

/*! Returns the sum-of-magnitudes of the elements in an n-D array (useful for L_1 norm) */
INLINE double ArraySumAbsnD(int     nvars,  /*!< number of elements at one array location, can be > 1 for systems of equations */
                            int     ndims,  /*!< number of dimensions */
                            int     *dim,   /*!< integer array of size in each dimension */
                            int     ghosts, /*!< number of ghost points in the array x */
                            int     *index, /*!< pre-allocated (by the calling function) integer array of size ndims */
                            double  *x      /*!< the array */
                           )
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

/*! Returns the sum-of-squares of the elements in an n-D array (useful for L_2 norm) */
INLINE double ArraySumSquarenD(int    nvars,  /*!< number of elements at one array location, can be > 1 for systems of equations */
                               int    ndims,  /*!< number of dimensions */
                               int    *dim,   /*!< integer array of size in each dimension */
                               int    ghosts, /*!< number of ghost points in the array x */
                               int    *index, /*!< pre-allocated (by the calling function) integer array of size ndims */
                               double *x      /*!< the array */
                              )
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
