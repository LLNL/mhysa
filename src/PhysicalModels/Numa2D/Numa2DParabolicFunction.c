#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d.h>
#include <mpivars.h>
#include <hypar.h>

/*
    These are not the actual viscous terms for compressible flows.
    Reference: Giraldo, Restelli, "A study of spectral element and discontinuous 
               Galerkin methods for the Navier–Stokes equations in nonhydrostatic 
               mesoscale atmospheric modeling: Equation sets and test cases", 
               Journal of Computational Physics, 227 (2008), pp. 3849--3877
*/

int Numa2DParabolicFunction(double *par,double *u,void *s,void *m,double t)
{
  HyPar           *solver   = (HyPar*) s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  Numa2D          *physics  = (Numa2D*) solver->physics;
  int             i,v,done;
  double          dxinv, dyinv;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int imax   = solver->dim_local[0];
  int jmax   = solver->dim_local[1];
  int *dim   = solver->dim_local;
  int size   = (imax+2*ghosts)*(jmax+2*ghosts)*_MODEL_NVARS_;

  _ArraySetValue_(par,size,0.0);
  if (physics->mu <= 0) return(0); /* inviscid flow */
  solver->count_par++;

  double        mu        = physics->mu;        

  /* allocate some arrays */
  double *Q, *QDeriv, *FViscous, *FDeriv; 
  Q         = (double*) calloc (size,sizeof(double)); /* primitive variables                */
  QDeriv    = (double*) calloc (size,sizeof(double)); /* derivative of primitive variables  */
  FViscous  = (double*) calloc (size,sizeof(double)); /* viscous flux                       */
  FDeriv    = (double*) calloc (size,sizeof(double)); /* derivative of viscous flux         */

  int index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];
  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); for (i=0; i<_MODEL_NDIMS_; i++) bounds[i] += 2*ghosts;
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double drho,uvel,vvel,dT,rho0,T0,P0,EP,ycoord;

    _GetCoordinate_             (_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->x,ycoord);
    physics->StandardAtmosphere (physics,ycoord,&EP,&P0,&rho0,&T0);
    _Numa2DGetFlowVars_         ( (u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);

    Q[_MODEL_NVARS_*p+0] = rho0 + drho;         /* density              */
    Q[_MODEL_NVARS_*p+1] = uvel;                /* x-velocity           */
    Q[_MODEL_NVARS_*p+2] = vvel;                /* y-velocity           */
    Q[_MODEL_NVARS_*p+3] = (dT+T0)/(drho+rho0); /* potential temperature*/

    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  /* Along X */
  _ArraySetValue_(QDeriv,size,0.0);
  IERR solver->FirstDerivativePar(QDeriv,Q,_XDIR_,1,solver,mpi);                    CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,dim,ghosts,mpi,QDeriv);CHECKERR(ierr);
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    _GetCoordinate_(_XDIR_,index[_XDIR_]-ghosts,dim,ghosts,solver->dxinv,dxinv);

    double rho, ux, vx, tx;
    rho = Q     [_MODEL_NVARS_*p+0];
    ux  = QDeriv[_MODEL_NVARS_*p+1] / dxinv;
    vx  = QDeriv[_MODEL_NVARS_*p+2] / dxinv;
    tx  = QDeriv[_MODEL_NVARS_*p+3] / dxinv;
    
    FViscous[_MODEL_NVARS_*p+0] = 0.0;
    FViscous[_MODEL_NVARS_*p+1] = mu * rho * ux;
    FViscous[_MODEL_NVARS_*p+2] = mu * rho * vx;
    FViscous[_MODEL_NVARS_*p+3] = mu * rho * tx;

    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }
  _ArraySetValue_(FDeriv,size,0.0);
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,-1,solver,mpi);             CHECKERR(ierr);
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
    for (v=0; v<_MODEL_NVARS_; v++) (par+p)[v] += (dxinv * (FDeriv+p)[v] );
    _ArrayIncrementIndex_(_MODEL_NDIMS_,dim,index,done);
  }

  /* Along Y */
  _ArraySetValue_(QDeriv,size,0.0);
  IERR solver->FirstDerivativePar(QDeriv,Q,_YDIR_,1,solver,mpi);                    CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,dim,ghosts,mpi,QDeriv);CHECKERR(ierr);
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->dxinv,dyinv);

    double rho, uy, vy, ty;
    rho = Q     [_MODEL_NVARS_*p+0];
    uy  = QDeriv[_MODEL_NVARS_*p+1] / dyinv;
    vy  = QDeriv[_MODEL_NVARS_*p+2] / dyinv;
    ty  = QDeriv[_MODEL_NVARS_*p+3] / dyinv;

    FViscous[_MODEL_NVARS_*p+0] = 0.0;
    FViscous[_MODEL_NVARS_*p+1] = mu * rho * uy;
    FViscous[_MODEL_NVARS_*p+2] = mu * rho * vy;
    FViscous[_MODEL_NVARS_*p+3] = mu * rho * ty;

    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }
  _ArraySetValue_(FDeriv,size,0.0);
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,-1,solver,mpi);             CHECKERR(ierr);
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
    for (v=0; v<_MODEL_NVARS_; v++) (par+p)[v] += (dyinv * (FDeriv+p)[v] );
    _ArrayIncrementIndex_(_MODEL_NDIMS_,dim,index,done);
  }

  free(Q);
  free(QDeriv);
  free(FViscous);
  free(FDeriv);

  return(0);
}
