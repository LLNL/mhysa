/*! @file NavierStokes3DParabolicFunction.c
    @author Debojyoti Ghosh
    @brief Compute the viscous terms for the 3D Navier Stokes equations
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Compute the viscous terms in the 3D Navier Stokes equations: this function computes
    the following:
    \f{equation}{
      \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ \tau_{zx} \\ u \tau_{xx} + v \tau_{yx} + w \tau_{zx} - q_x \end{array}\right]
      + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{zy} \\ u \tau_{xy} + v \tau_{yy} + w \tau_{zy} - q_y \end{array}\right]
      + \frac {\partial} {\partial z} \left[\begin{array}{c} 0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\ u \tau_{xz} + v \tau_{yz} + w \tau_{zz} - q_z \end{array}\right]
    \f}
    where 
    \f{align}{
      \tau_{xx} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(2\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} - \frac{\partial w}{\partial z}\right),\\
      \tau_{xy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{xz} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right),\\
      \tau_{yx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yy} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} +2\frac{\partial v}{\partial y} - \frac{\partial w}{\partial z}\right),\\
      \tau_{yz} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial w}{\partial y} + \frac{\partial v}{\partial z}\right),\\
      \tau_{zx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right),\\
      \tau_{zy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}\right),\\
      \tau_{zz} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} + 2\frac{\partial v}{\partial y}\right),\\
      q_x &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial x}, \\
      q_y &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial y}, \\
      q_z &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial z}
    \f}
    and the temperature is \f$T = \gamma p/\rho\f$. \f$Re\f$ and \f$Pr\f$ are the Reynolds and Prandtl numbers, respectively. Note that this function
    computes the entire parabolic term, and thus bypasses HyPar's parabolic function calculation interfaces. NavierStokes3DInitialize() assigns this
    function to #HyPar::ParabolicFunction.
    \n\n
    Reference:
    + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
      Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
      number is based on speed of sound, instead of the freestream velocity).
*/
int NavierStokes3DParabolicFunction(
                                      double  *par, /*!< Array to hold the computed viscous terms */
                                      double  *u,   /*!< Solution vector array */
                                      void    *s,   /*!< Solver object of type #HyPar */
                                      void    *m,   /*!< MPI object of type #MPIVariables */
                                      double  t     /*!< Current simulation time */
                                   )
{
  HyPar           *solver   = (HyPar*) s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  NavierStokes3D  *physics  = (NavierStokes3D*) solver->physics;
  int             i,j,k,v;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int imax   = solver->dim_local[0];
  int jmax   = solver->dim_local[1];
  int kmax   = solver->dim_local[2];
  int *dim   = solver->dim_local;
  int size   = solver->npoints_local_wghosts;

  _ArraySetValue_(par,size*_MODEL_NVARS_,0.0);
  if (physics->Re <= 0) return(0); /* inviscid flow */
  solver->count_par++;

  static double two_third    = 2.0/3.0;
  double        inv_gamma_m1 = 1.0 / (physics->gamma-1.0);
  double        inv_Re       = 1.0 / physics->Re;
  double        inv_Pr       = 1.0 / physics->Pr;

  double *Q; /* primitive variables */
  Q = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double energy,pressure;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
        _NavierStokes3DGetFlowVar_( (u+p),
                                    Q[p+0],
                                    Q[p+1],
                                    Q[p+2],
                                    Q[p+3],
                                    energy,
                                    pressure,
                                    physics);
        Q[p+4] = physics->gamma*pressure/Q[p+0]; /* temperature */
      }
    }
  }

  double *QDerivX = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));
  double *QDerivY = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));
  double *QDerivZ = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));

  IERR solver->FirstDerivativePar(QDerivX,Q,_XDIR_,1,solver,mpi); CHECKERR(ierr);
  IERR solver->FirstDerivativePar(QDerivY,Q,_YDIR_,1,solver,mpi); CHECKERR(ierr);
  IERR solver->FirstDerivativePar(QDerivZ,Q,_ZDIR_,1,solver,mpi); CHECKERR(ierr);

  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->dim_local,
                                 solver->ghosts,mpi,QDerivX); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->dim_local,
                                 solver->ghosts,mpi,QDerivY); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->dim_local,
                                 solver->ghosts,mpi,QDerivY); CHECKERR(ierr);

  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dxinv, dyinv, dzinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
        _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
        _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
        _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->dxinv,dzinv);
        _ArrayScale1D_((QDerivX+p),dxinv,_MODEL_NVARS_);
        _ArrayScale1D_((QDerivY+p),dyinv,_MODEL_NVARS_);
        _ArrayScale1D_((QDerivZ+p),dzinv,_MODEL_NVARS_);
      }
    }
  }

  double *FViscous = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));
  double *FDeriv   = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));

  /* Along X */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;

        double uvel, vvel, wvel, T, Tx, 
               ux, uy, uz, vx, vy, wx, wz;
        uvel = (Q+p)[1];
        vvel = (Q+p)[2];
        wvel = (Q+p)[3];
        T    = (Q+p)[4];
        Tx   = (QDerivX+p)[4];
        ux   = (QDerivX+p)[1];
        vx   = (QDerivX+p)[2];
        wx   = (QDerivX+p)[3];
        uy   = (QDerivY+p)[1];
        vy   = (QDerivY+p)[2];
        uz   = (QDerivZ+p)[1];
        wz   = (QDerivZ+p)[3];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);

        double tau_xx, tau_xy, tau_xz, qx;
        tau_xx = two_third * (mu*inv_Re) * (2*ux - vy - wz);
        tau_xy = (mu*inv_Re) * (uy + vx);
        tau_xz = (mu*inv_Re) * (uz + wx);
        qx     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tx;

        (FViscous+p)[0] = 0.0;
        (FViscous+p)[1] = tau_xx;
        (FViscous+p)[2] = tau_xy;
        (FViscous+p)[3] = tau_xz;
        (FViscous+p)[4] = uvel*tau_xx + vvel*tau_xy + wvel*tau_xz + qx;
      }
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,-1,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dxinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
        _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
        for (v=0; v<_MODEL_NVARS_; v++) (par+p)[v] += (dxinv * (FDeriv+p)[v] );
      }
    }
  }

  /* Along Y */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;

        double uvel, vvel, wvel, T, Ty, 
               ux, uy, vx, vy, vz, wy, wz;
        uvel = (Q+p)[1];
        vvel = (Q+p)[2];
        wvel = (Q+p)[3];
        T    = (Q+p)[4];
        Ty   = (QDerivY+p)[4];
        ux   = (QDerivX+p)[1];
        vx   = (QDerivX+p)[2];
        uy   = (QDerivY+p)[1];
        vy   = (QDerivY+p)[2];
        wy   = (QDerivY+p)[3];
        vz   = (QDerivZ+p)[2];
        wz   = (QDerivZ+p)[3];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);

        double tau_yx, tau_yy, tau_yz, qy;
        tau_yx = (mu*inv_Re) * (uy + vx);
        tau_yy = two_third * (mu*inv_Re) * (-ux + 2*vy - wz);
        tau_yz = (mu*inv_Re) * (vz + wy);
        qy     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Ty;

        (FViscous+p)[0] = 0.0;
        (FViscous+p)[1] = tau_yx;
        (FViscous+p)[2] = tau_yy;
        (FViscous+p)[3] = tau_yz;
        (FViscous+p)[4] = uvel*tau_yx + vvel*tau_yy + wvel*tau_yz + qy;
      }
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,-1,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dyinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
        _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
        for (v=0; v<_MODEL_NVARS_; v++) (par+p)[v] += (dyinv * (FDeriv+p)[v] );
      }
    }
  }

  /* Along Z */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;

        double uvel, vvel, wvel, T, Tz, 
               ux, uz, vy, vz, wx, wy, wz;
        uvel = (Q+p)[1];
        vvel = (Q+p)[2];
        wvel = (Q+p)[3];
        T    = (Q+p)[4];
        Tz   = (QDerivZ+p)[4];
        ux   = (QDerivX+p)[1];
        wx   = (QDerivX+p)[3];
        vy   = (QDerivY+p)[2];
        wy   = (QDerivY+p)[3];
        uz   = (QDerivZ+p)[1];
        vz   = (QDerivZ+p)[2];
        wz   = (QDerivZ+p)[3];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T,0.76);

        double tau_zx, tau_zy, tau_zz, qz;
        tau_zx = (mu*inv_Re) * (uz + wx);
        tau_zy = (mu*inv_Re) * (vz + wy);
        tau_zz = two_third * (mu*inv_Re) * (-ux - vy + 2*wz);
        qz     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tz;

        (FViscous+p)[0] = 0.0;
        (FViscous+p)[1] = tau_zx;
        (FViscous+p)[2] = tau_zy;
        (FViscous+p)[3] = tau_zz;
        (FViscous+p)[4] = uvel*tau_zx + vvel*tau_zy + wvel*tau_zz + qz;
      }
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_ZDIR_,-1,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dzinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
        _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->dxinv,dzinv);
        for (v=0; v<_MODEL_NVARS_; v++) (par+p)[v] += (dzinv * (FDeriv+p)[v] );
      }
    }
  }

  free(Q);
  free(QDerivX);
  free(QDerivY);
  free(QDerivZ);
  free(FViscous);
  free(FDeriv);

  if (solver->flag_ib) _ArrayBlockMultiply_(par,solver->iblank,size,_MODEL_NVARS_);
  return(0);
}
