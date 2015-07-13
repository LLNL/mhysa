/*! @file BCSWSlipWall.c
    @author Debojyoti Ghosh
    @brief Slip-wall boundary conditions for shallow water equations
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/shallowwater1d.h>
#include <physicalmodels/shallowwater2d.h>

/*! Applies the slip-wall boundary condition: This is specific to the one and two dimenstional
    shallow water equations (#ShallowWater1D, #ShallowWater2D).
    It is used for simulating inviscid walls or symmetric boundaries. The height,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary 
    face is equal to the specified wall velocity.
*/
int BCSWSlipWallU(
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

  if (ndims == 1) {

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
        
        /* flow variables in the interior */
        double h, uvel;
        double h_gpt, uvel_gpt;
        _ShallowWater1DGetFlowVar_((phi+nvars*p2),h,uvel);
        /* set the ghost point values */
        h_gpt = h;
        uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;

        phi[nvars*p1+0] = h_gpt;
        phi[nvars*p1+1] = h_gpt * uvel_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  } else if (ndims == 2) {

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
        
        /* flow variables in the interior */
        double h, uvel, vvel;
        double h_gpt, uvel_gpt, vvel_gpt;
        _ShallowWater2DGetFlowVar_((phi+nvars*p2),h,uvel,vvel);
        /* set the ghost point values */
        h_gpt = h;
        if (dim == _XDIR_) {
          uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
          vvel_gpt = vvel;
        } else if (dim == _YDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
        } else {
          uvel_gpt = 0.0;
          vvel_gpt = 0.0;
        }

        phi[nvars*p1+0] = h_gpt;
        phi[nvars*p1+1] = h_gpt * uvel_gpt;
        phi[nvars*p1+2] = h_gpt * vvel_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}

/*! Applies the slip-wall boundary condition to the delta-solution: This is specific to the one
    and two dimensional shallow water equations (#ShallowWater1D, #ShallowWater2D).
    It is used for simulating inviscid walls or symmetric boundaries. The height,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary 
    face is equal to the specified wall velocity.
    \n\n
    The above treatment is applied on the delta-solution added to the reference solution.
*/
int BCSWSlipWallDU(
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
  int v;

  if (ndims == 1) {

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

        /* flow in the interior is phi + phi_ref (since phi is DU) */
        double phi_total[nvars]; 
        for (v=0; v<nvars; v++) phi_total[v] = phi[nvars*p2+v]+phi_ref[nvars*p2+v];
        
        /* flow variables in the interior */
        double h , uvel;
        double h0, uvel0;
        _ShallowWater1DGetFlowVar_(phi_total,h,uvel);
        _ShallowWater1DGetFlowVar_((phi_ref+nvars*p2),h0,uvel0);
        /* set the ghost point values */
        double h_gpt, uvel_gpt;
        double h0_gpt, uvel0_gpt;
        /* ghost point values of the total flow variables */
        h_gpt = h;
        uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
        /* ghost point value of the reference flow variables */
        h0_gpt = h0;
        uvel0_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel0;

        phi[nvars*p1+0] = h_gpt             - h0_gpt;
        phi[nvars*p1+1] = h_gpt * uvel_gpt  - h0_gpt * uvel0_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  } else if (ndims == 2) {

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

        /* flow in the interior is phi + phi_ref (since phi is DU) */
        double phi_total[nvars]; 
        for (v=0; v<nvars; v++) phi_total[v] = phi[nvars*p2+v]+phi_ref[nvars*p2+v];
        
        /* flow variables in the interior */
        double h , uvel , vvel;
        double h0, uvel0, vvel0;
        _ShallowWater2DGetFlowVar_(phi_total,h,uvel,vvel);
        _ShallowWater2DGetFlowVar_((phi_ref+nvars*p2),h0,uvel0,vvel0);
        /* set the ghost point values */
        double h_gpt, uvel_gpt, vvel_gpt;
        double h0_gpt, uvel0_gpt, vvel0_gpt;
        /* ghost point values of the total flow variables */
        h_gpt = h;
        if (dim == _XDIR_) {
          uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
          vvel_gpt = vvel;
        } else if (dim == _YDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
        } else {
          uvel_gpt = 0.0;
          vvel_gpt = 0.0;
        }
        /* ghost point value of the reference flow variables */
        h0_gpt = h0;
        if (dim == _XDIR_) {
          uvel0_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel0;
          vvel0_gpt = vvel0;
        } else if (dim == _YDIR_) {
          uvel0_gpt = uvel0;
          vvel0_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel0;
        } else {
          uvel0_gpt = 0.0;
          vvel0_gpt = 0.0;
        }

        phi[nvars*p1+0] = h_gpt             - h0_gpt;
        phi[nvars*p1+1] = h_gpt * uvel_gpt  - h0_gpt * uvel0_gpt;
        phi[nvars*p1+2] = h_gpt * vvel_gpt  - h0_gpt * vvel0_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}
