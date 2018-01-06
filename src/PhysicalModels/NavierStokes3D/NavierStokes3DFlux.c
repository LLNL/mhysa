/*! @file NavierStokes3DFlux.c
    @author Debojyoti Ghosh
    @brief Functions to compute the hyperbolic flux for 3D Navier-Stokes equations
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*!
  Compute the hyperbolic flux function for the 3D Navier-Stokes equations:
  Note: the flux function needs to be computed at the ghost points as well.
*/
int NavierStokes3DFlux(
                        double  *f, /*!< Array to hold the computed flux vector (same layout as u) */
                        double  *u, /*!< Array with the solution vector */
                        int     dir,/*!< Spatial dimension (x, y, or z) for which to compute the flux */
                        void    *s, /*!< Solver object of type #HyPar */
                        double  t   /*!< Current simulation time */
                      )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;
  int               i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ns      = param->n_species;
  int nv      = param->n_vibeng;

  static int index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho_s[ns], rho_t, vx, vy, vz, E, E_v[nv], P, T;
    _NavierStokes3DGetFlowVar_((u+nvars*p),rho_s,rho_t,vx,vy,vz,E,E_v,P,T,param);
    _NavierStokes3DSetFlux_((f+nvars*p),rho_s,rho_t,vx,vy,vz,E,E_v,P,param,dir);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  
  return(0);
}

