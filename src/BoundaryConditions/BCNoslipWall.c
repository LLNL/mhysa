#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/navierstokes2d.h>
#include <physicalmodels/navierstokes3d.h>

/* 
 * No-Slip Wall BC - specific to NavierStokes2D/NavierStokes3D
 * Used for viscous walls
 *
 * boundary->var is irrelevant, it acts on on the components
 * so no need to specify it for each component, just specify
 * it once with an arbitrary value for boundary->var 
*/

int BCNoslipWallU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (ndims == 2) {

    /* create a fake physics object */
    NavierStokes2D physics; 
    double gamma;
    gamma = physics.gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

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
        double rho, uvel, vvel, energy, pressure;
        double rho_gpt, uvel_gpt, vvel_gpt, energy_gpt, pressure_gpt;
        _NavierStokes2DGetFlowVar_((phi+nvars*p2),rho,uvel,vvel,energy,pressure,(&physics));
        /* set the ghost point values */
        rho_gpt = rho;
        pressure_gpt = pressure;
        uvel_gpt = 2.0*boundary->FlowVelocity[0] - uvel;
        vvel_gpt = 2.0*boundary->FlowVelocity[1] - vvel;
        energy_gpt = inv_gamma_m1*pressure_gpt
                    + 0.5 * rho_gpt * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt);

        phi[nvars*p1+0] = rho_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt;
        phi[nvars*p1+3] = energy_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  } else if (ndims == 3) {

    /* create a fake physics object */
    NavierStokes3D physics;
    double gamma;
    gamma = physics.gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

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
        double rho, uvel, vvel, wvel, energy, pressure;
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho,uvel,vvel,wvel,energy,pressure,(&physics));
        /* set the ghost point values */
        rho_gpt = rho;
        pressure_gpt = pressure;
        uvel_gpt = 2.0*boundary->FlowVelocity[0] - uvel;
        vvel_gpt = 2.0*boundary->FlowVelocity[1] - vvel;
        wvel_gpt = 2.0*boundary->FlowVelocity[2] - wvel;
        energy_gpt = inv_gamma_m1*pressure_gpt
                    + 0.5 * rho_gpt 
                    * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);

        phi[nvars*p1+0] = rho_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt;
        phi[nvars*p1+3] = rho_gpt * wvel_gpt;
        phi[nvars*p1+4] = energy_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}

int BCNoslipWallDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

  if (ndims == 2) {

    /* create a fake physics object */
    NavierStokes2D physics; 
    double gamma; 
    gamma = physics.gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

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
        double rho , uvel , vvel , energy , pressure ;
        double rho0, uvel0, vvel0, energy0, pressure0;
        _NavierStokes2DGetFlowVar_(phi_total,rho,uvel,vvel,energy,pressure,(&physics));
        _NavierStokes2DGetFlowVar_((phi_ref+nvars*p2),rho0,uvel0,vvel0,energy0,pressure0,(&physics));
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, energy_gpt, pressure_gpt;
        double rho0_gpt, uvel0_gpt, vvel0_gpt, energy0_gpt, pressure0_gpt;
        /* ghost point values of the total flow variables */
        rho_gpt = rho;
        pressure_gpt = pressure;
        uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
        vvel_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
        energy_gpt = inv_gamma_m1*pressure_gpt 
                    + 0.5 * rho_gpt * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt);
        /* ghost point value of the reference flow variables */
        rho0_gpt = rho0;
        pressure0_gpt = pressure0;
        uvel0_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel0;
        vvel0_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel0;
        energy0_gpt = inv_gamma_m1*pressure0_gpt 
                    + 0.5 * rho0_gpt * (uvel0_gpt*uvel0_gpt + vvel0_gpt*vvel0_gpt);

        phi[nvars*p1+0] = rho_gpt             - rho0_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt  - rho0_gpt * uvel0_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt  - rho0_gpt * vvel0_gpt;
        phi[nvars*p1+3] = energy_gpt          - energy0_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  } else if (ndims == 3) {

    /* create a fake physics object */
    NavierStokes3D physics; 
    double gamma; 
    gamma = physics.gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

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
        double rho, uvel, vvel, wvel, energy, pressure;
        double rho0, uvel0, vvel0, wvel0, energy0, pressure0;
        _NavierStokes3DGetFlowVar_(phi_total,rho,uvel,vvel,wvel,energy,pressure,(&physics));
        _NavierStokes3DGetFlowVar_((phi_ref+nvars*p2),rho0,uvel0,vvel0,wvel0,energy0,pressure0,(&physics));
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        double rho0_gpt, uvel0_gpt, vvel0_gpt, wvel0_gpt, energy0_gpt, pressure0_gpt;
        /* setting the ghost point values for the total flow variables */
        rho_gpt = rho;
        pressure_gpt = pressure;
        uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
        vvel_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
        wvel_gpt = 2.0*boundary->FlowVelocity[_ZDIR_] - wvel;
        energy_gpt = inv_gamma_m1*pressure_gpt 
                    + 0.5 * rho_gpt 
                    * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);
        /* setting the ghost point values for the reference flow variables */
        rho0_gpt = rho0;
        pressure0_gpt = pressure0;
        uvel0_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel0;
        vvel0_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel0;
        wvel0_gpt = 2.0*boundary->FlowVelocity[_ZDIR_] - wvel0;
        energy0_gpt = inv_gamma_m1*pressure0_gpt 
                    + 0.5 * rho0_gpt 
                    * (uvel0_gpt*uvel0_gpt + vvel0_gpt*vvel0_gpt + wvel0_gpt*wvel0_gpt);

        phi[nvars*p1+0] = rho_gpt            - rho0_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt - rho0_gpt * uvel0_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt - rho0_gpt * vvel0_gpt;
        phi[nvars*p1+3] = rho_gpt * wvel_gpt - rho0_gpt * wvel0_gpt;
        phi[nvars*p1+4] = energy_gpt         - energy0_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}
