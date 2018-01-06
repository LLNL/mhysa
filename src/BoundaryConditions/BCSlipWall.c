/*! @file BCSlipWall.c
    @author Debojyoti Ghosh
    @brief Slip-wall boundary conditions
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/euler1d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Applies the slip-wall boundary condition: This is specific to the 1D Euler (#Euler1D) and 
    3D Navier-Stokes system (#NavierStokes3D).
    It is used for simulating inviscid walls or symmetric boundaries. The pressure, density,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary 
    face is equal to the specified wall velocity.
*/
int BCSlipWallU(
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

  int i;

  if (ndims == 1) {

    Euler1D *physics = (Euler1D*) (*(boundary->physics)); 
    int ns = physics->n_species;
    int nv = physics->n_vibeng;
    double gamma = physics->gamma; 
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
        
        double rho_s[ns], rho_t, uvel, E, E_v[nv], P, T;
        _Euler1DGetFlowVar_((phi+nvars*p2),rho_s,rho_t,uvel,E,E_v,P,T,physics);

        double rho_s_g[ns], rho_t_g, uvel_g, E_g, E_v_g[nv], P_g;
        for (i = 0; i < ns; i++) rho_s_g[i] = rho_s[i];
        rho_t_g = rho_t;
        uvel_g = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
        P_g = P;
        E_g = inv_gamma_m1*P_g/rho_t_g + 0.5*uvel_g*uvel_g;
        for (i = 0; i < nv; i++) E_v_g[i] = E_v[i];
        _Euler1DSetFlowVar_((phi+nvars*p1),rho_s_g,rho_t_g,uvel_g,E_g,E_v_g,P_g,physics);

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  } else if (ndims == 3) {

    NavierStokes3D *physics = (NavierStokes3D*) (*(boundary->physics)); 
    int ns = physics->n_species;
    int nv = physics->n_vibeng;
    double gamma = physics->gamma; 
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
        
        double rho_s[ns], rho_t, uvel, vvel, wvel, E, E_v[nv], P, T;
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho_s,rho_t,uvel,vvel,wvel,E,E_v,P,T,physics);

        double rho_s_g[ns], rho_t_g, uvel_g, vvel_g, wvel_g, E_g, E_v_g[nv], P_g;
        for (i = 0; i < ns; i++) rho_s_g[i] = rho_s[i];
        rho_t_g = rho_t;
        if (dim == _XDIR_) {
          uvel_g = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
          vvel_g = vvel;
          wvel_g = wvel;
        } else if (dim == _YDIR_) {
          uvel_g = uvel;
          vvel_g = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
          wvel_g = wvel;
        } else if (dim == _ZDIR_) {
          uvel_g = uvel;
          vvel_g = vvel;
          wvel_g = 2.0*boundary->FlowVelocity[_ZDIR_] - wvel;
        } else {
          uvel_g = 0.0;
          vvel_g = 0.0;
          wvel_g = 0.0;
        }
        P_g = P;
        E_g = inv_gamma_m1*P_g/rho_t_g + 0.5*(uvel_g*uvel_g+vvel_g*vvel_g+wvel_g*wvel_g);
        for (i = 0; i < nv; i++) E_v_g[i] = E_v[i];
        _NavierStokes3DSetFlowVar_((phi+nvars*p1),rho_s_g,rho_t_g,uvel_g,vvel_g,wvel_g,E_g,E_v_g,P_g,physics);

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}
