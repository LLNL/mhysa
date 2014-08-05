/* second derivative scheme definitions */
#define _SECOND_ORDER_CENTRAL_  "2"
#define _FOURTH_ORDER_CENTRAL_  "4"

/* One-dimensional second derivative functions: Functions to calculate the
 * finite-difference approximation to the second derivative along a given 
 * dimension at the cell-centers / grid points.
 *
 * Arguments:-
 *  
 *  Df        double*     array containing the approximation to the second derivative
 *                        (1D array representing an n-D solution)
 *
 *  f         double*     array containing the function whose second derivative is to
 *                        be computed.
 *                        (1D array representing an n-D solution)
 *
 *  dir       int         dimension (x/y/z/...) along which to compute the second derivative
 *
 *  s         void*       pointer to an object providing the solver context. The object must
 *                        contain at least the following:
 *                        ndims     int       number of dimensions
 *                        nvars     int       number of elements at each grid location
 *                        dim       int*      array of grid size along each dimension 
 *                                            (integer array of size ndims)
 *                        ghosts    int       number of ghost points
 *
 *  m         void*       Pointer to an object providing the MPI context. Current not being 
 *                        used. This will be required for compact finite-difference methods.
 *
 *
 * Notes:
 *
 *  + The second derivative is computed at the ghost points too. Thus, biased schemes are used
 *    at and near the boundaries.
*/

/* Second derivative functions */
int SecondDerivativeSecondOrderCentral (double*,double*,int,void*,void*);
int SecondDerivativeFourthOrderCentral (double*,double*,int,void*,void*);
