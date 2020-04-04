/*! @file NavierStokes3DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 3D Navier Stokes equations.
*/
#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*! Compute the Roe-averaged state for the 3D Navier Stokes equations. This function 
    just calls the macro #_NavierStokes3DRoeAverage_ and is not used by any 
    functions within the 3D Navier Stokes module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes3DRoeAverage(
                              double  *uavg, /*!< The computed Roe-averaged state */
                              double  *uL,   /*!< Left state (conserved variables)*/
                              double  *uR,   /*!< Right state (conserved variables)*/
                              void    *p     /*!< Object of type #NavierStokes3D with physics-related variables */
                            )
{
  NavierStokes3D *param  = (NavierStokes3D*) p;
  _NavierStokes3DRoeAverage_(uavg,uL,uR,param);
  return(0);
}

/*! Compute the pressure from the conserved solution on a grid */
int NavierStokes3DComputePressure(  double*             P, /*!< Array to hold the computed pressure (same layout as u) */
                                    const double* const u, /*!< Array with the solution vector */
                                    void*               s  /*!< Solver object of type #HyPar */
                                  )
{
  HyPar          *solver = (HyPar*)   s;
  NavierStokes3D *param  = (NavierStokes3D*) solver->physics;
  int            i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ns      = param->n_species;
  int nv      = param->n_vibeng;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int idx; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,idx);
    double rho_s[ns], rho_t, vx, vy, vz, E, E_v[nv], tmp, T;
    _NavierStokes3DGetFlowVar_((u+nvars*idx),rho_s,rho_t,vx,vy,vz,E,E_v,tmp,T,param);
    P[idx] = tmp;
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  
  return(0);
}

/*! Compute the temperature from the conserved solution on a grid */
int NavierStokes3DComputeTemperature( double*             T, /*!< Array to hold the computed pressure (same layout as u) */
                                      const double* const u, /*!< Array with the solution vector */
                                      void*               s  /*!< Solver object of type #HyPar */
                                    )
{
  HyPar          *solver = (HyPar*)   s;
  NavierStokes3D *param  = (NavierStokes3D*) solver->physics;
  int            i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ns      = param->n_species;
  int nv      = param->n_vibeng;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int idx; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,idx);
    double rho_s[ns], rho_t, vx, vy, vz, E, E_v[nv], P, tmp;
    _NavierStokes3DGetFlowVar_((u+nvars*idx),rho_s,rho_t,vx,vy,vz,E,E_v,P,tmp,param);
    T[idx] = tmp;
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  
  return(0);
}
