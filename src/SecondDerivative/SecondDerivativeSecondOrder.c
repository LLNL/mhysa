/*! @file SecondDerivativeSecondOrder.c
    @brief 2nd order discretization of the second derivative.
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <secondderivative.h>
#include <mpivars.h>
#include <hypar.h>

#undef  _MINIMUM_GHOSTS_
/*! \def _MINIMUM_GHOSTS_
 * Minimum number of ghost points required.
*/
#define _MINIMUM_GHOSTS_ 1

/*! Computes the second-order finite-difference approximation to the second derivative (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial^2 f\right)_i = f_{i-1} - 2f_i + f_{i+1}
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative.

    \b Notes:
    + The second derivative is computed at the grid points or the cell centers.
    + Though the array D2f includes ghost points, the second derivative is \b not computed at these 
      locations. Thus, array elements corresponding to the ghost points contain undefined values.
    + \a D2f and \a f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.
*/
int SecondDerivativeSecondOrderCentral(
                                        double  *D2f, /*!< Array to hold the computed second derivative (with ghost points)
                                                           (same size and layout as f) */
                                        double  *f,   /*!< Array containing the grid point function values whose first
                                                           derivative is to be computed (with ghost points) */
                                        int     dir,  /*!< The spatial dimension along which the derivative is computed */
                                        void    *s,   /*!< Solver object of type #HyPar */
                                        void    *m    /*!< MPI object of type #MPIVariables */
                                      )
{
  HyPar         *solver = (HyPar*) s;
  int           i, v;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  if ((!D2f) || (!f)) {
    fprintf(stderr, "Error in SecondDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }
  if (ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in SecondDerivativeSecondOrder(): insufficient number of ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    for (i = 0; i < dim[dir]; i++) {
      int qL, qC, qR;
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC);
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  D2f[qC*nvars+v] = f[qL*nvars+v]-2*f[qC*nvars+v]+f[qR*nvars+v];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }
  
  return(0);
}
