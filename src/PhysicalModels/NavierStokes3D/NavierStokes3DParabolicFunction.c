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
      \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ \tau_{xx} \\ \tau_{yx} \\ \tau_{zx} \\ u \tau_{xx} + v \tau_{yx} + w \tau_{zx} - q_x \\ 0 \\ \vdots \\ 0 \end{array}\right]
      + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{zy} \\ u \tau_{xy} + v \tau_{yy} + w \tau_{zy} - q_y \\ 0 \\ \vdots \\ 0 \end{array}\right]
      + \frac {\partial} {\partial z} \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\ u \tau_{xz} + v \tau_{yz} + w \tau_{zz} - q_z \\ 0 \\ \vdots \\ 0 \end{array}\right]
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

  int nvars = solver->nvars;
  int ns    = physics->n_species;
  int nv    = physics->n_vibeng;

  int ghosts = solver->ghosts;
  int imax   = solver->dim_local[0];
  int jmax   = solver->dim_local[1];
  int kmax   = solver->dim_local[2];
  int *dim   = solver->dim_local;
  int size   = solver->npoints_local_wghosts;

  _ArraySetValue_(par,size*nvars,0.0);
  if (physics->Re <= 0) return(0); /* inviscid flow */
  solver->count_par++;

  static double two_third    = 2.0/3.0;
  double        inv_gamma_m1 = 1.0 / (physics->gamma-1.0);
  double        inv_Re       = 1.0 / physics->Re;
  double        inv_Pr       = 1.0 / physics->Pr;

  double *Q; /* primitive variables */
  Q = (double*) calloc (size*nvars,sizeof(double));

  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {

        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;

        double rho_s[ns],rho_t, uvel, vvel, wvel, E, E_v[nv], P, T;
        _NavierStokes3DGetFlowVar_( (u+p),rho_s,rho_t,uvel,vvel,wvel,E,E_v,P,T,physics);
        for ( v = 0; v < ns; v++) Q[p+v] = rho_s[v];
        Q[p+ns]   = uvel;
        Q[p+ns+1] = vvel;
        Q[p+ns+2] = wvel;
        Q[p+ns+3] = physics->gamma*P/rho_t;
        for ( v = 0; v < nv; v++) Q[p+ns+4+v] = E_v[v];

      }
    }
  }

  double *QDerivX = (double*) calloc (size*nvars,sizeof(double));
  double *QDerivY = (double*) calloc (size*nvars,sizeof(double));
  double *QDerivZ = (double*) calloc (size*nvars,sizeof(double));

  IERR solver->FirstDerivativePar(QDerivX,Q,_XDIR_,1,solver,mpi); CHECKERR(ierr);
  IERR solver->FirstDerivativePar(QDerivY,Q,_YDIR_,1,solver,mpi); CHECKERR(ierr);
  IERR solver->FirstDerivativePar(QDerivZ,Q,_ZDIR_,1,solver,mpi); CHECKERR(ierr);

  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,nvars,solver->dim_local,
                                 solver->ghosts,mpi,QDerivX); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,nvars,solver->dim_local,
                                 solver->ghosts,mpi,QDerivY); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,nvars,solver->dim_local,
                                 solver->ghosts,mpi,QDerivY); CHECKERR(ierr);

  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {

        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dxinv, dyinv, dzinv;

        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;
        _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
        _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
        _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->dxinv,dzinv);

        _ArrayScale1D_((QDerivX+p),dxinv,nvars);
        _ArrayScale1D_((QDerivY+p),dyinv,nvars);
        _ArrayScale1D_((QDerivZ+p),dzinv,nvars);

      }
    }
  }

  double *FViscous = (double*) calloc (size*nvars,sizeof(double));
  double *FDeriv   = (double*) calloc (size*nvars,sizeof(double));

  /* Along X */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;

        double uvel, vvel, wvel, T, Tx, 
               ux, uy, uz, vx, vy, wx, wz;
        uvel = (Q+p)[ns];
        vvel = (Q+p)[ns+1];
        wvel = (Q+p)[ns+2];
        T    = (Q+p)[ns+3];
        Tx   = (QDerivX+p)[ns+3];
        ux   = (QDerivX+p)[ns];
        vx   = (QDerivX+p)[ns+1];
        wx   = (QDerivX+p)[ns+2];
        uy   = (QDerivY+p)[ns];
        vy   = (QDerivY+p)[ns+1];
        uz   = (QDerivZ+p)[ns];
        wz   = (QDerivZ+p)[ns+2];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);

        double tau_xx, tau_xy, tau_xz, qx;
        tau_xx = two_third * (mu*inv_Re) * (2*ux - vy - wz);
        tau_xy = (mu*inv_Re) * (uy + vx);
        tau_xz = (mu*inv_Re) * (uz + wx);
        qx     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tx;

        for (v = 0; v < ns; v++) (FViscous+p)[v] = 0.0;
        (FViscous+p)[ns]    = tau_xx;
        (FViscous+p)[ns+1]  = tau_xy;
        (FViscous+p)[ns+2]  = tau_xz;
        (FViscous+p)[ns+3]  = uvel*tau_xx + vvel*tau_xy + wvel*tau_xz + qx;
        for (v = 0; v < nv; v++) (FViscous+p)[ns+4+v]  = 0.0;
      }
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,-1,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dxinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;
        _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
        for (v=0; v<nvars; v++) (par+p)[v] += (dxinv * (FDeriv+p)[v] );
      }
    }
  }

  /* Along Y */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;

        double uvel, vvel, wvel, T, Ty, 
               ux, uy, vx, vy, vz, wy, wz;
        uvel = (Q+p)[ns];
        vvel = (Q+p)[ns+1];
        wvel = (Q+p)[ns+2];
        T    = (Q+p)[ns+3];
        Ty   = (QDerivY+p)[ns+3];
        ux   = (QDerivX+p)[ns];
        vx   = (QDerivX+p)[ns+1];
        uy   = (QDerivY+p)[ns];
        vy   = (QDerivY+p)[ns+1];
        wy   = (QDerivY+p)[ns+2];
        vz   = (QDerivZ+p)[ns+1];
        wz   = (QDerivZ+p)[ns+2];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);

        double tau_yx, tau_yy, tau_yz, qy;
        tau_yx = (mu*inv_Re) * (uy + vx);
        tau_yy = two_third * (mu*inv_Re) * (-ux + 2*vy - wz);
        tau_yz = (mu*inv_Re) * (vz + wy);
        qy     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Ty;

        for (v = 0; v < ns; v++) (FViscous+p)[v] = 0.0;
        (FViscous+p)[ns]    = tau_yx;
        (FViscous+p)[ns+1]  = tau_yy;
        (FViscous+p)[ns+2]  = tau_yz;
        (FViscous+p)[ns+3]  = uvel*tau_yx + vvel*tau_yy + wvel*tau_yz + qy;
        for (v = 0; v < nv; v++) (FViscous+p)[ns+4+v]  = 0.0;
      }
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,-1,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dyinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;
        _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
        for (v=0; v<nvars; v++) (par+p)[v] += (dyinv * (FDeriv+p)[v] );
      }
    }
  }

  /* Along Z */
  for (i=-ghosts; i<(imax+ghosts); i++) {
    for (j=-ghosts; j<(jmax+ghosts); j++) {
      for (k=-ghosts; k<(kmax+ghosts); k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;

        double uvel, vvel, wvel, T, Tz, 
               ux, uz, vy, vz, wx, wy, wz;
        uvel = (Q+p)[ns];
        vvel = (Q+p)[ns+1];
        wvel = (Q+p)[ns+2];
        T    = (Q+p)[ns+3];
        Tz   = (QDerivZ+p)[ns+3];
        ux   = (QDerivX+p)[ns];
        wx   = (QDerivX+p)[ns+2];
        vy   = (QDerivY+p)[ns+1];
        wy   = (QDerivY+p)[ns+2];
        uz   = (QDerivZ+p)[ns];
        vz   = (QDerivZ+p)[ns+1];
        wz   = (QDerivZ+p)[ns+2];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T,0.76);

        double tau_zx, tau_zy, tau_zz, qz;
        tau_zx = (mu*inv_Re) * (uz + wx);
        tau_zy = (mu*inv_Re) * (vz + wy);
        tau_zz = two_third * (mu*inv_Re) * (-ux - vy + 2*wz);
        qz     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tz;

        for (v = 0; v < ns; v++) (FViscous+p)[v] = 0.0;
        (FViscous+p)[ns]    = tau_zx;
        (FViscous+p)[ns+1]  = tau_zy;
        (FViscous+p)[ns+2]  = tau_zz;
        (FViscous+p)[ns+3]  = uvel*tau_zx + vvel*tau_zy + wvel*tau_zz + qz;
        for (v = 0; v < nv; v++) (FViscous+p)[ns+4+v]  = 0.0;
      }
    }
  }
  IERR solver->FirstDerivativePar(FDeriv,FViscous,_ZDIR_,-1,solver,mpi); CHECKERR(ierr);
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p,index[3]; index[0]=i; index[1]=j; index[2]=k;
        double dzinv;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= nvars;
        _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->dxinv,dzinv);
        for (v=0; v<nvars; v++) (par+p)[v] += (dzinv * (FDeriv+p)[v] );
      }
    }
  }

  free(Q);
  free(QDerivX);
  free(QDerivY);
  free(QDerivZ);
  free(FViscous);
  free(FDeriv);

  if (solver->flag_ib) _ArrayBlockMultiply_(par,solver->iblank,size,nvars);
  return(0);
}
