/*! @file BCSponge.c
    @author Debojyoti Ghosh
    @brief Sponge boundary condition
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

/*! Applies the sponge boundary condition: This function computes the source term required to apply a
    sponge boundary condition that gradually relaxes the solution to a specified state. This boundary
    condition is different from other boundary conditions in the sense that it is applied at interior
    grid points (but within the defined sponge zone).
    \n\n
    The source term for the sponge is computed as:
    \f{align}{
      {\bf S}_i &= \sigma_i \left( {\bf U}_i - {\bf U}_{\rm ref} \right),\\
      \sigma_i &= \frac {x_i - x_{\rm start}} {x_{\rm end} - x_{\rm start}}
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the sponge, \f${\bf U}_{\rm ref}\f$
    is the specified state to which the solution is relaxed, and \f$x_i\f$, \f$x_{\rm start}\f$, and 
    \f$x_{\rm end}\f$ are the spatial coordinates of the grid point, sponge start, and sponge end, 
    respectively, along the spatial dimension of the sponge.
*/
int BCSpongeSource(
                    void    *b,     /*!< Boundary object of type #DomainBoundary */
                    int     ndims,  /*!< Number of spatial dimensions */
                    int     nvars,  /*!< Number of variables/DoFs per grid point */
                    int     ghosts, /*!< Number of ghost points */
                    int     *size,  /*!< Integer array with the number of grid points in each spatial dimension */
                    double  *grid,  /*!< 1D array with the spatial coordinates of the grid points, one dimension after the other */
                    double  *u,     /*!< Solution */
                    double  *source /*!< Source term to which the sponge term is added */
                  )
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            dim       = boundary->dim;
  int            face      = boundary->face;
  double         *uref     = boundary->SpongeValue;
  double         *xmin     = boundary->xmin;
  double         *xmax     = boundary->xmax;
  int            v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int i = indexb[dim] + boundary->is[dim];
      double x, xstart, xend;
      _GetCoordinate_(dim,i,size,ghosts,grid,x);
      xstart = xmin[dim];
      xend   = xmax[dim];
      /* calculate sigma */
      double sigma;
      if (face > 0) sigma = (x - xstart) / (xend - xstart);
      else          sigma = (x - xend  ) / (xstart - xend);
      /* add to the source term */
      int p; _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p);
      for (v=0; v<nvars; v++) source[nvars*p+v] -= (sigma * (u[nvars*p+v]-uref[v]));
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

/*! Dummy function to ensure consistency with the overall boundary condition 
    implementation. The actual sponge boundary condition is implemented by
    BCSpongeSource()
*/
int BCSpongeUDummy(
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
  return(0);
}
