#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

/*
  Refer: Computational Fluid Mechanics and Heat Transfer
         by Tannehill, Anderson and Pletcher
         Chapter 5, Section 5.1.7 for the non-dimensional
         form of the NS equations.
*/

int NavierStokes2DParabolicFunction(double *par,double *u,void *s,void *m,double t)
{
  HyPar           *solver   = (HyPar*) s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  NavierStokes2D  *physics  = (NavierStokes2D*) solver->physics;
  int             i,j,v,d;
  _DECLARE_IERR_;


  int ghosts = solver->ghosts;
  int imax   = solver->dim_local[0];
  int jmax   = solver->dim_local[1];
  int *dim   = solver->dim_local;
  int nvars  = solver->nvars;
  int ndims  = solver->ndims;
  int size   = (imax+2*ghosts)*(jmax+2*ghosts)*nvars;

  static double two_third = 2.0/3.0;

  _ArraySetValue_(par,size,0.0);
  if (physics->Re <= 0) return(0); /* inviscid flow */

  double *Q; /* primitive variables */
  Q = (double*) calloc (size,sizeof(double));
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      int p,index[2]; index[0]=i; index[1]=j;
      double energy,pressure;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;
      _NavierStokes2DGetFlowVar_( (u+p),Q[p+0],Q[p+1],Q[p+2],energy,
                                  pressure,physics);
      Q[p+3] = physics->gamma*(physics->Minf*physics->Minf)*pressure/Q[p+0]; /* temperature */
    }
  }

  double *QDerivX = (double*) calloc (size,sizeof(double));
  double *QDerivY = (double*) calloc (size,sizeof(double));
  IERR solver->FirstDerivativePar(QDerivX,Q,_XDIR_,solver,mpi); CHECKERR(ierr);
  IERR solver->FirstDerivativePar(QDerivY,Q,_YDIR_,solver,mpi); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,dim,
                               solver->ghosts,mpi,QDerivX);     CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,dim,
                               solver->ghosts,mpi,QDerivY);     CHECKERR(ierr);
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      int p,index[2]; index[0]=i; index[1]=j;
      double dxinv, dyinv;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;
      _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
      _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
      _ArrayScale1D_((QDerivX+p),dxinv,nvars);
      _ArrayScale1D_((QDerivY+p),dyinv,nvars);
    }
  }

  double *FViscous = (double*) calloc (size,sizeof(double));
  double *FDeriv   = (double*) calloc (size,sizeof(double));

  /* Along X */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      int p,index[2]; index[0]=i; index[1]=j;;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;

      double uvel, vvel, T, Tx, ux, uy, vx, vy;
      uvel = (Q+p)[1];
      vvel = (Q+p)[2];
      T    = (Q+p)[3];
      Tx   = (QDerivX+p)[3];
      ux   = (QDerivX+p)[1];
      vx   = (QDerivX+p)[2];
      uy   = (QDerivY+p)[1];
      vy   = (QDerivY+p)[2];

      /* calculate viscosity coeff based on Sutherland's law */
      /* double mu = physics->C1 * raiseto(T,1.5) / (T + physics->C2); */
      double mu = 1.0; 

      double tau_xx, tau_xy, qx;
      tau_xx = two_third * (mu/physics->Re) * (2*ux - vy);
      tau_xy = (mu/physics->Re) * (uy + vx);
      qx     = ( (mu/physics->Re) * (1.0/(physics->gamma-1.0)) * (1.0/physics->Pr) ) * (1.0/(physics->Minf*physics->Minf)) * Tx;

      (FViscous+p)[0] = 0.0;
      (FViscous+p)[1] = tau_xx;
      (FViscous+p)[2] = tau_xy;
      (FViscous+p)[3] = uvel*tau_xx + vvel*tau_xy + qx;
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      int p,index[2]; index[0]=i; index[1]=j;
      double dxinv;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;
      _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
      for (v=0; v<nvars; v++) (par+p)[v] += (dxinv * (FDeriv+p)[v] );
    }
  }

  /* Along Y */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      int p,index[2]; index[0]=i; index[1]=j;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;

      double uvel, vvel, T, Ty, ux, uy, vx, vy;
      uvel = (Q+p)[1];
      vvel = (Q+p)[2];
      T    = (Q+p)[3];
      Ty   = (QDerivY+p)[3];
      ux   = (QDerivX+p)[1];
      vx   = (QDerivX+p)[2];
      uy   = (QDerivY+p)[1];
      vy   = (QDerivY+p)[2];

      /* calculate viscosity coeff based on Sutherland's law */
      /* double mu = physics->C1 * raiseto(T,1.5) / (T + physics->C2); */
      double mu = 1.0; 

      double tau_yx, tau_yy, qy;
      tau_yx = (mu/physics->Re) * (uy + vx);
      tau_yy = two_third * (mu/physics->Re) * (-ux + 2*vy);
      qy     = ( (mu/physics->Re) * (1.0/(physics->gamma-1.0)) * (1.0/physics->Pr) ) * (1.0/(physics->Minf*physics->Minf)) * Ty;

      (FViscous+p)[0] = 0.0;
      (FViscous+p)[1] = tau_yx;
      (FViscous+p)[2] = tau_yy;
      (FViscous+p)[3] = uvel*tau_yx + vvel*tau_yy + qy;
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      int p,index[2]; index[0]=i; index[1]=j;
      double dyinv;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;
      _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
      for (v=0; v<nvars; v++) (par+p)[v] += (dyinv * (FDeriv+p)[v] );
    }
  }

  free(Q);
  free(QDerivX);
  free(QDerivY);
  free(FViscous);
  free(FDeriv);

  return(0);
}
