/* first derivative scheme definitions */
#define _SECOND_ORDER_CENTRAL_  "2"
#define _FOURTH_ORDER_CENTRAL_  "4"

/* One-dimensional first derivative functions: Functions to calculate the
 * finite-difference approximation to the first derivative along a given 
 * dimension at the cell-centers / grid points.
 *
 * Arguments:-
 *  
 *  Df        double*     array containing the approximation to the first derivative
 *                        (1D array representing an n-D solution)
 *
 *  f         double*     array containing the function whose first derivative is to
 *                        be computed.
 *                        (1D array representing an n-D solution)
 *
 *  dir       int         dimension (x/y/z/...) along which to compute the first derivative
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
 *  + The first derivative is computed at the ghost points too. Thus, biased schemes are used
 *    at and near the boundaries.
*/

/* First derivative functions */
int FirstDerivativeSecondOrderCentral (double*,double*,int,void*,void*);
int FirstDerivativeFourthOrderCentral (double*,double*,int,void*,void*);
