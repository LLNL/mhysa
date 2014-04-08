#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes3d.h>

/* 
 * Slip Wall BC - specific to Navier-Stokes 3D
 *
 * boundary->var is irrelevant, it acts on on the components
 * so no need to specify it for each component, just specify
 * it once with an arbitrary value for boundary->var 
*/

int BCSlipWall(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (ndims == 2) {

    /* create a fake physics object */
    Euler2D physics; physics.gamma = 1.4;

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
        else if (face == -1) indexi[dim] = size[dim]-indexb[dim];
        else return(1);
        _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
        
        /* flow variables in the interior */
        double rho, uvel, vvel, energy, pressure;
        _Euler2DGetFlowVar_((phi+nvars*p2),rho,uvel,vvel,energy,pressure,(&physics));
        /* set the ghost point values */
        phi[nvars*p1+0] = rho;
        if (dim == _XDIR_) {
          phi[nvars*p1+1] = rho * (-uvel);
          phi[nvars*p1+2] = rho * ( vvel);
        } else if (dim == _YDIR_) {
          phi[nvars*p1+1] = rho * ( uvel);
          phi[nvars*p1+2] = rho * (-vvel);
        }
        phi[nvars*p1+3] = energy;
        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  } else if (ndims == 3) {

    /* create a fake physics object */
    NavierStokes3D physics; physics.gamma = 1.4;

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
        else if (face == -1) indexi[dim] = size[dim]-indexb[dim];
        else return(1);
        _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
        
        /* flow variables in the interior */
        double rho, uvel, vvel, wvel, energy, pressure;
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho,uvel,vvel,wvel,energy,pressure,(&physics));
        /* set the ghost point values */
        phi[nvars*p1+0] = rho;
        if (dim == _XDIR_) {
          phi[nvars*p1+1] = rho * (-uvel);
          phi[nvars*p1+2] = rho * ( vvel);
          phi[nvars*p1+3] = rho * ( wvel);
        } else if (dim == _YDIR_) {
          phi[nvars*p1+1] = rho * ( uvel);
          phi[nvars*p1+2] = rho * (-vvel);
          phi[nvars*p1+3] = rho * ( wvel);
        } else if (dim == _ZDIR_) {
          phi[nvars*p1+1] = rho * ( uvel);
          phi[nvars*p1+2] = rho * ( vvel);
          phi[nvars*p1+3] = rho * (-wvel);
        }
        phi[nvars*p1+4] = energy;
        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}

