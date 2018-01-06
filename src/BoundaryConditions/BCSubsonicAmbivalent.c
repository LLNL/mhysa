/*! @file BCSubsonicAmbivalent.c
    @author Debojyoti Ghosh
    @brief Subsonic "ambivalent" boundary conditions (specific to Euler/Navier-Stokes systems)
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/navierstokes3d.h>

/*! Applies the subsonic "ambivalent" boundary condition: 
    The velocity at the boundary face is extrapolated from the interior and 
    its dot product with the boundary normal (pointing into the domain) is
    computed. 
    + If it is positive, subsonic inflow boundary conditions are applied: the 
      density and velocity at the physical boundary ghost points are specified, 
      while the pressure is extrapolated from the interior of the domain. 
    + If it is negative, subsonic outflow boundary conditions are applied: the 
      pressure at the physical boundary ghost points is specified, while the 
      density and velocity are extrapolated from the interior. 
        
    This boundary condition is specific to the 3D Navier-Stokes systems 
    (#NavierStokes3D).
*/
int BCSubsonicAmbivalentU(
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
  int k;

  if (ndims == 3) {

    NavierStokes3D *physics = (NavierStokes3D*) (*(boundary->physics));
    double gamma = physics->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);
    int ns = physics->n_species;
    int nv = physics->n_vibeng;

    /* boundary normal (pointing into the domain) */
    double nx, ny, nz;
    if (dim == 0) {
      nx = 1.0;
      ny = 0.0;
      nz = 0.0;
    } else if (dim == 1) {
      nx = 0.0;
      ny = 1.0;
      nz = 0.0;
    } else if (dim == 2) {
      nx = 0.0;
      ny = 0.0;
      nz = 1.0;
    }
    nx *= (double) face;
    ny *= (double) face;
    nz *= (double) face;

    if (boundary->on_this_proc) {

      int bounds[ndims], indexb[ndims], indexi[ndims], indexj[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);

      int done = 0;
      while (!done) {

        int p1, p2;
        double rho_s[ns], rho_t, uvel, vvel, wvel, E, E_v[nv], pressure, T;

        /* compute boundary face velocity  - 2nd order */
        _ArrayCopy1D_(indexb,indexi,ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
        _ArrayCopy1D_(indexi,indexj,ndims);
        if (face ==  1) {
          indexi[dim] = 0;
          indexj[dim] = indexi[dim] + 1;
        } else if (face == -1) {
          indexi[dim] = size[dim]-1;
          indexj[dim] = indexi[dim] - 1;
        }
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexj,ghosts,p2);

        double uvel1, uvel2, uvelb,
               vvel1, vvel2, vvelb,
               wvel1, wvel2, wvelb;
        _NavierStokes3DGetFlowVar_((phi+nvars*p1),rho_s,rho_t,uvel1,vvel1,wvel1,E,E_v,pressure,T,physics);
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho_s,rho_t,uvel2,vvel2,wvel2,E,E_v,pressure,T,physics);
        uvelb = 1.5*uvel1 - 0.5*uvel2;
        vvelb = 1.5*vvel1 - 0.5*vvel2;
        wvelb = 1.5*wvel1 - 0.5*wvel2;
        double vel_normal = uvelb*nx + vvelb*ny + wvelb*nz;

        _ArrayCopy1D_(indexb,indexi,ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
        if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
        
        /* flow variables in the interior */
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho_s,rho_t,uvel,vvel,wvel,E,E_v,pressure,T,physics);

        /* set the ghost point values */
        double rho_s_gpt[ns], rho_t_gpt, uvel_gpt, vvel_gpt, wvel_gpt, E_gpt, E_v_gpt[nv], pressure_gpt;
        if (vel_normal > 0) {
          /* inflow */
          for (k = 0; k < ns; k++)  rho_s_gpt[k] = boundary->FlowDensity[k];
          _NavierStokes3DTotalDensity_(rho_t_gpt, rho_s_gpt, ns);
          pressure_gpt = pressure;
          uvel_gpt = boundary->FlowVelocity[0];
          vvel_gpt = boundary->FlowVelocity[1];
          wvel_gpt = boundary->FlowVelocity[2];
        } else {
          /* outflow */
          for (k = 0; k < ns; k++) rho_s_gpt[k] = rho_s[k];
          rho_t_gpt = rho_t;
          pressure_gpt = boundary->FlowPressure;
          uvel_gpt = uvel;
          vvel_gpt = vvel;
          wvel_gpt = wvel;
        }
        E_gpt = inv_gamma_m1*pressure_gpt/rho_t_gpt
                + 0.5 * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);
        for (k = 0; k < nv; k++) E_v_gpt[k] = E_v[k];

        _NavierStokes3DSetFlowVar_( (phi+nvars*p1),rho_s_gpt,rho_t_gpt,
                                    uvel_gpt,vvel_gpt,wvel_gpt,
                                    E_gpt,E_v_gpt,pressure_gpt,physics );

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  
  return(0);
}
