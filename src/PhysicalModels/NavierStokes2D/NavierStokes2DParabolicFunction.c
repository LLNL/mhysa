/*! @file NavierStokes2DParabolicFunction.c
    @author Debojyoti Ghosh
    @brief Compute the viscous terms for the 2D Navier Stokes equations
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Compute the viscous terms in the 2D Navier Stokes equations: this function computes
    the following:
    \f{equation}{
      \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ u \tau_{xx} + v \tau_{yx} - q_x \end{array}\right]
      + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ u \tau_{xy} + v \tau_{yy} - q_y \end{array}\right]
    \f}
    where 
    \f{align}{
      \tau_{xx} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(2\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}\right),\\
      \tau_{xy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yy} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} +2\frac{\partial v}{\partial y}\right),\\
      q_x &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial x}, \\
      q_y &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial y}
    \f}
    and the temperature is \f$T = \gamma p/\rho\f$. \f$Re\f$ and \f$Pr\f$ are the Reynolds and Prandtl numbers, respectively. Note that this function
    computes the entire parabolic term, and thus bypasses HyPar's parabolic function calculation interfaces. NavierStokes2DInitialize() assigns this
    function to #HyPar::ParabolicFunction.
    \n\n
    Reference:
    + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
      Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
      number is based on speed of sound, instead of the freestream velocity).
*/
int NavierStokes2DParabolicFunction(
                                      double  *par, /*!< Array to hold the computed viscous terms */
                                      double  *u,   /*!< Solution vector array */
                                      void    *s,   /*!< Solver object of type #HyPar */
                                      void    *m,   /*!< MPI object of type #MPIVariables */
                                      double  t     /*!< Current simulation time */
                                   )
{
  HyPar           *solver   = (HyPar*) s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  NavierStokes2D  *physics  = (NavierStokes2D*) solver->physics;
  int             i,j,v;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int imax   = solver->dim_local[0];
  int jmax   = solver->dim_local[1];
  int *dim   = solver->dim_local;
  int nvars  = solver->nvars;
  int ndims  = solver->ndims;
  int size   = (imax+2*ghosts)*(jmax+2*ghosts)*nvars;

  _ArraySetValue_(par,size,0.0);
  if (physics->Re <= 0) return(0); /* inviscid flow */
  solver->count_par++;

  static double two_third = 2.0/3.0;
  double        inv_gamma_m1 = 1.0 / (physics->gamma-1.0);
  double        inv_Re       = 1.0 / physics->Re;
  double        inv_Pr       = 1.0 / physics->Pr;

  double *Q; /* primitive variables */
  Q = (double*) calloc (size,sizeof(double));
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      int p,index[2]; index[0]=i; index[1]=j;
      double energy,pressure;
      _ArrayIndex1D_(ndims,dim,index,ghosts,p); p *= nvars;
      _NavierStokes2DGetFlowVar_( (u+p),Q[p+0],Q[p+1],Q[p+2],energy,
                                  pressure,physics);
      Q[p+3] = physics->gamma*pressure/Q[p+0]; /* temperature */
    }
  }

  double *QDerivX = (double*) calloc (size,sizeof(double));
  double *QDerivY = (double*) calloc (size,sizeof(double));

  IERR solver->FirstDerivativePar(QDerivX,Q,_XDIR_,1,solver,mpi); CHECKERR(ierr);
  IERR solver->FirstDerivativePar(QDerivY,Q,_YDIR_,1,solver,mpi); CHECKERR(ierr);

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
      double mu = raiseto(T, 0.76);

      double tau_xx, tau_xy, qx;
      tau_xx = two_third * (mu*inv_Re) * (2*ux - vy);
      tau_xy = (mu*inv_Re) * (uy + vx);
      qx     = ( (mu*inv_Re) * inv_gamma_m1 * inv_Pr ) * Tx;

      (FViscous+p)[0] = 0.0;
      (FViscous+p)[1] = tau_xx;
      (FViscous+p)[2] = tau_xy;
      (FViscous+p)[3] = uvel*tau_xx + vvel*tau_xy + qx;
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,-1,solver,mpi); CHECKERR(ierr);
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
      double mu = raiseto(T, 0.76);

      double tau_yx, tau_yy, qy;
      tau_yx = (mu*inv_Re) * (uy + vx);
      tau_yy = two_third * (mu*inv_Re) * (-ux + 2*vy);
      qy     = ( (mu*inv_Re) * inv_gamma_m1 * inv_Pr ) * Ty;

      (FViscous+p)[0] = 0.0;
      (FViscous+p)[1] = tau_yx;
      (FViscous+p)[2] = tau_yy;
      (FViscous+p)[3] = uvel*tau_yx + vvel*tau_yy + qy;
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,-1,solver,mpi); CHECKERR(ierr);
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
