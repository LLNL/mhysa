/*! @file BCDirichlet.c
    @author Debojyoti Ghosh
    @brief Dirichlet boundary conditions
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

/*! Applies (steady) Dirichlet boundary conditions for the solution: the ghost points at the physical
    boundaries are set to specified values */
int BCDirichletU(
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

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,size  ,indexb,boundary->is,ghosts,p);
      _ArrayCopy1D_((boundary->DirichletValue),(phi+nvars*p),nvars);
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

/*! Applies (steady) Dirichlet boundary conditions for the delta-solution: the ghost points at the physical
    boundaries are set to zero (since the change in the steady Dirichlet boundary state is zero */
int BCDirichletDU(
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

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,size  ,indexb,boundary->is,ghosts,p);
      _ArraySetValue_((phi+nvars*p),nvars,0.0);
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}
