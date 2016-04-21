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
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ \rho u w \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ \rho v w \\ (e+p)v \end{array}\right], \\
    dir = z, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho w \\ \rho u w \\ \rho v w \\ \rho w^2 + p \\ (e+p)w \end{array}\right]
  \f}
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
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, vx, vy, vz, e, P;
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vx,vy,vz,e,P,param);
    _NavierStokes3DSetFlux_((f+_MODEL_NVARS_*p),rho,vx,vy,vz,e,P,dir);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  
  return(0);
}

/*! Compute the stiff flux, given the solution vector. The stiff flux is the component of
    the total flux that represents the acoustic modes (see #_NavierStokes3DSetStiffFlux_).
    Here, the linearized approximation to the stiff flux is computed as:
    \f{equation}{
      {\bf f}\left({\bf u}\right) = A_f{\bf u}
    \f}
    where \f$A_f = A_f\left({\bf u}_{ref}\right)\f$ is the fast Jacobian (#NavierStokes3D::fast_jac)
    evaluated for the solution at the beginning of each time step (\f${\bf u}_{ref}\f$ is 
    #NavierStokes3D::solution). This is done in NavierStokes3DPreStep().\n\n
  Note: the flux function needs to be computed at the ghost points as well.
*/
int NavierStokes3DStiffFlux(
                              double  *f, /*!< Array to hold the computed flux vector (same layout as u) */
                              double  *u, /*!< Array with the solution vector */
                              int     dir,/*!< Spatial dimension (x,y, or z) for which to compute the flux */
                              void    *s, /*!< Solver object of type #HyPar */
                              double  t   /*!< Current simulation time */
                           )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double *Af = param->fast_jac+(_MODEL_NDIMS_*p+dir)*JacSize;
    MatVecMult5(_MODEL_NVARS_,(f+_MODEL_NVARS_*p),Af,(u+_MODEL_NVARS_*p)); 
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}

/*! Compute the non-stiff flux, given the solution vector. The non-stiff flux is the component of
    the total flux that represents the entropy modes (see #_NavierStokes3DSetNonStiffFlux_).
    Here, the linearized approximation to the non-stiff flux is computed as:
    \f{equation}{
      {\bf f}\left({\bf u}\right) - A_f{\bf u}
    \f}
    where \f${\bf f}\left({\bf u}\right)\f$ is the total flux computed in NavierStokes3DFlux(), 
    and \f$A_f{\bf u}\f$ is the linearized stiff flux computed in NavierStokes3DStiffFlux().\n\n
  Note: the flux function needs to be computed at the ghost points as well.
*/
int NavierStokes3DNonStiffFlux(
                              double  *f, /*!< Array to hold the computed flux vector (same layout as u) */
                              double  *u, /*!< Array with the solution vector */
                              int     dir,/*!< Spatial dimension (x,y, or z) for which to compute the flux */
                              void    *s, /*!< Solver object of type #HyPar */
                              double  t   /*!< Current simulation time */
                           )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];
  static double     ftot[_MODEL_NVARS_], fstiff[_MODEL_NVARS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    /* compute total flux */
    double rho, vx, vy, vz, e, P;
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vx,vy,vz,e,P,param);
    _NavierStokes3DSetFlux_(ftot,rho,vx,vy,vz,e,P,dir);
    /* compute stiff stuff */
    double *Af = param->fast_jac+(_MODEL_NDIMS_*p+dir)*JacSize;
    MatVecMult5(_MODEL_NVARS_,fstiff,Af,(u+_MODEL_NVARS_*p)); 
    /* subtract stiff flux from total flux */
    _ArraySubtract1D_((f+_MODEL_NVARS_*p),ftot,fstiff,_MODEL_NVARS_);
    /* Done */
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}
