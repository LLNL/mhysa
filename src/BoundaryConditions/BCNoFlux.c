/*! @file BCNoFlux.c
    @author Debojyoti Ghosh
    @brief No-flux boundary condition (specific to #Numa2D and #Numa3D).
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/numa2d.h>
#include <physicalmodels/numa3d.h>

/*! Applies the no-flux boundary conditions: This boundary condition is specific 
    to the NUMA 2D/3D (#Numa2D, #Numa3D). Used for simulating inviscid walls or 
    symmetry boundaries. It's equivalent to the slip-wall BC of the Euler/Navier-
    Stokes system.\n\n
    The density, potential temperature, and tangential velocity are extrapolated, 
    while the normal velocity at the ghost point is set to the negative of that in 
    the interior (to enforce zero-normal velocity at the boundary face).
*/
int BCNoFluxU(
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
    _ArraySubtract1D_ (bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_   (indexb,ndims,0);
    int done = 0;
    while (!done) {
      int p1, p2;
      _ArrayCopy1D_ (indexb,indexi,ndims);
      _ArrayAdd1D_  (indexi,indexi,boundary->is,ndims);
      if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
      else return(1);
      _ArrayIndex1DWO_  (ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_    (ndims,size,indexi,ghosts,p2);
      
      if (nvars == 4) {
        phi[nvars*p1+0] = phi[nvars*p2+0];
        phi[nvars*p1+1] = (dim == _XDIR_ ? -phi[nvars*p2+1] : phi[nvars*p2+1] );
        phi[nvars*p1+2] = (dim == _YDIR_ ? -phi[nvars*p2+2] : phi[nvars*p2+2] );
        phi[nvars*p1+3] = phi[nvars*p2+3];
      } else if (nvars == 5) {
        phi[nvars*p1+0] = phi[nvars*p2+0];
        phi[nvars*p1+1] = (dim == _XDIR_ ? -phi[nvars*p2+1] : phi[nvars*p2+1] );
        phi[nvars*p1+2] = (dim == _YDIR_ ? -phi[nvars*p2+2] : phi[nvars*p2+2] );
        phi[nvars*p1+3] = (dim == _ZDIR_ ? -phi[nvars*p2+3] : phi[nvars*p2+3] );
        phi[nvars*p1+4] = phi[nvars*p2+4];
      }

      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

/*! Applies the no-flux boundary conditions to the delta-solution by setting the physical boundary ghost 
    point delta-solution values to zero.
*/
int BCNoFluxDU(
                void    *b,       /*!< Boundary object of type #DomainBoundary */
                void    *m,       /*!< MPI object of type #MPIVariables */
                int     ndims,    /*!< Number of spatial dimensions */
                int     nvars,    /*!< Number of variables/DoFs per grid point */
                int     *size,    /*!< Integer array with the number of grid points in each spatial dimension */
                int     ghosts,   /*!< Number of ghost points */
                double  *phi,     /*!< The solution array on which to apply the boundary condition -
                                       Note that this is a delta-solution \f$\Delta {\bf U}\f$.*/
                double  *phi_ref, /*!< Reference solution */
                double  waqt      /*!< Current solution time */
              )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    _ArraySubtract1D_ (bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_   (indexb,ndims,0);
    int done = 0;
    while (!done) {
      int p1, p2;
      _ArrayCopy1D_ (indexb,indexi,ndims);
      _ArrayAdd1D_  (indexi,indexi,boundary->is,ndims);
      if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
      else return(1);
      _ArrayIndex1DWO_  (ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_    (ndims,size,indexi,ghosts,p2);
      _ArraySetValue_   ((phi+nvars*p1),nvars,0.0);
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }

  return(0);
}
