/*! @file BCReflect.c
    @author Debojyoti Ghosh
    @brief Reflection boundary condition
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

/*! Applies the reflection boundary condition: The values at the physical 
    boundary ghost points are set to the negative of the interior values
    adjacent to the boundary.
*/
int BCReflectU(
                void    *b,     /*!< Boundary object of type #DomainBoundary */
                void    *m,     /*!< MPI object of type #MPIVariables */
                int     ndims,  /*!< Number of spatial dimensions */
                int     nvars,  /*!< Number of variables/DoFs per grid point */
                int     *size,  /*!< Integer array with the number of grid points in each spatial dimension */
                int     ghosts, /*!< Number of ghost points */
                double  *phi,   /*!< The solution array on which to apply the boundary condition */
                double  waqt    /*!< Current solution time */
              )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0);
    int done = 0;
    while (!done) {
      int p1, p2;
      _ArrayCopy1D_(indexb,indexi,ndims);
      _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
      if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
      else return(1);
      _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
      _ArrayScaleCopy1D_((phi+nvars*p2),(-1.0),(phi+nvars*p1),nvars);
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}
