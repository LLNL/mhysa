/*! @file BCPeriodic.c
    @author Debojyoti Ghosh
    @brief Periodic boundary conditions
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>

/*! Applies periodic boundary conditions: Implemented by copying the solution 
    from the other end of the domain into the physical boundary ghost points.
    \n\n
    **Note**: This function only acts if the the number of processors is 1 along
    the spatial dimension this boundary corresponds to. If there are more than 1
    processors along this dimension, periodicity is handled by MPIExchangeBoundariesnD()
    to minimize communication.
*/
int BCPeriodicU(
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
  MPIVariables   *mpi      = (MPIVariables*)   m;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if ((boundary->on_this_proc) && (mpi->iproc[dim] == 1)) {
    int bounds[ndims], index1[ndims], index2[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(index1,ndims,0);
    _ArraySetValue_(index2,ndims,0);
    int done = 0;
    while (!done) {
      int p1 = 0, p2 = 0;
      _ArrayCopy1D_(index1,index2,ndims);
      if (face == 1) {
        index2[dim] = index1[dim] + size[dim]-ghosts;
        _ArrayIndex1DWO_(ndims,size,index1,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,index2,ghosts,p2);
      } else if (face == -1) {
        _ArrayIndex1DWO_(ndims,size,index1,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,index1,ghosts,p2);
      }
      _ArrayCopy1D_((phi+nvars*p2),(phi+nvars*p1),nvars);
      _ArrayIncrementIndex_(ndims,bounds,index1,done);
    }
  }
  return(0);
}
