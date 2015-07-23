/*! @file NavierStokes2DFlux.c
    @author Debojyoti Ghosh
    @brief Functions to compute the hyperbolic flux for 2D Navier-Stokes equations
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*!
  Compute the hyperbolic flux function for the 2D Navier-Stokes equations:
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ (e+p)v \end{array}\right]
  \f}
*/
int NavierStokes2DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar           *solver = (HyPar*)   s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             *dim    = solver->dim_local;
  int             ghosts  = solver->ghosts;
  static int      index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double rho, vx, vy, e, P;
    _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);
    _NavierStokes2DSetFlux_((f+_MODEL_NVARS_*p),rho,vx,vy,e,P,dir);
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}

int NavierStokes2DStiffFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes2D    *param  = (NavierStokes2D*) solver->physics;
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
    double *Af = param->fast_jac+(2*p+dir)*JacSize;
    MatVecMult4(_MODEL_NVARS_,(f+_MODEL_NVARS_*p),Af,(u+_MODEL_NVARS_*p)); 
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}
