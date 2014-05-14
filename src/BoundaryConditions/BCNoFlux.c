#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/numa3d.h>

/* 
 * No-Flux BC - specific to the NUMA 3D system
 * Used for inviscid walls or symmetry boundaries
 * (It's equivalent to the slip-wall BC of the Euler/
 * Navier-Stokes system, but the solution vector 
 * is (drho,u,v,w,dtheta) where theta is the temperature
 * potential.
 *
 * boundary->var is irrelevant, it acts on on the components
 * so no need to specify it for each component, just specify
 * it once with an arbitrary value for boundary->var 
*/

int BCNoFluxU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double waqt)
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
      
      /* flow variables in the interior */
      double drho, uvel, vvel, wvel, dT;
      double drho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, dT_gpt;
      _Numa3DGetFlowVars_((phi+nvars*p2),drho,uvel,vvel,wvel,dT,1.0);
      /* set the ghost point values */
      drho_gpt  = drho;
      dT_gpt    = dT;
      if (dim == _XDIR_) {
        uvel_gpt = -uvel;
        vvel_gpt =  vvel;
        wvel_gpt =  wvel;
      } else if (dim == _YDIR_) {
        uvel_gpt =  uvel;
        vvel_gpt = -vvel;
        wvel_gpt =  wvel;
      } else if (dim == _ZDIR_) {
        uvel_gpt =  uvel;
        vvel_gpt =  vvel;
        wvel_gpt = -wvel;
      }

      phi[nvars*p1+0] = drho_gpt;
      phi[nvars*p1+1] = uvel_gpt;
      phi[nvars*p1+2] = vvel_gpt;
      phi[nvars*p1+3] = wvel_gpt;
      phi[nvars*p1+4] = dT_gpt;

      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

int BCNoFluxDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

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
      
      /* flow variables in the interior */
      double drho, uvel, vvel, wvel, dT;
      double drho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, dT_gpt;
      _Numa3DGetFlowVars_((phi+nvars*p2),drho,uvel,vvel,wvel,dT,1.0);
      /* setting the ghost point values for the total flow variables */
      drho_gpt  = drho;
      dT_gpt    = dT;
      if (dim == _XDIR_) {
        uvel_gpt = -uvel;
        vvel_gpt =  vvel;
        wvel_gpt =  wvel;
      } else if (dim == _YDIR_) {
        uvel_gpt =  uvel;
        vvel_gpt = -vvel;
        wvel_gpt =  wvel;
      } else if (dim == _ZDIR_) {
        uvel_gpt =  uvel;
        vvel_gpt =  vvel;
        wvel_gpt = -wvel;
      }
      phi[nvars*p1+0] = drho_gpt;
      phi[nvars*p1+1] = uvel_gpt;
      phi[nvars*p1+2] = vvel_gpt;
      phi[nvars*p1+3] = wvel_gpt;
      phi[nvars*p1+4] = dT_gpt;

      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }

  return(0);
}
