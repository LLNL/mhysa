#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

inline int Euler2DGetFlowVar        (double*,double*,double*,double*,double*,double*,void*);
inline int Euler2DRoeAverage        (double*,double*,double*,void*);
inline int Euler2DEigenvalues       (double*,double**,void*,int);
inline int Euler2DLeftEigenvectors  (double*,double**,void*,int);
inline int Euler2DRightEigenvectors (double*,double**,void*,int);

int Euler2DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      ierr    = 0,done,k,v;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *index_inter  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  int *bounds_inter = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  double **R, **D, **L, **DL,**modA;
  R     = (double**) calloc (nvars,sizeof(double*));
  D     = (double**) calloc (nvars,sizeof(double*));
  L     = (double**) calloc (nvars,sizeof(double*));
  DL    = (double**) calloc (nvars,sizeof(double*));
  modA  = (double**) calloc (nvars,sizeof(double*));
  for (k = 0; k < nvars; k++){
    R[k]    = (double*) calloc (nvars,sizeof(double));
    D[k]    = (double*) calloc (nvars,sizeof(double));
    L[k]    = (double*) calloc (nvars,sizeof(double));
    DL[k]   = (double*) calloc (nvars,sizeof(double));
    modA[k] = (double*) calloc (nvars,sizeof(double));
  }

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double udiff[5], uavg[5],udiss[5];

      /* Roe's upwinding scheme */

      for (v=0; v<nvars; v++) udiff[v] = 0.5 * (uR[nvars*p+v] - uL[nvars*p+v]);

      ierr = Euler2DRoeAverage(&uavg[0],&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      ierr = Euler2DEigenvalues       (&uavg[0],D,param,dir); CHECKERR(ierr);
      ierr = Euler2DLeftEigenvectors  (&uavg[0],L,param,dir); CHECKERR(ierr);
      ierr = Euler2DRightEigenvectors (&uavg[0],R,param,dir); CHECKERR(ierr);

      for (k=0; k<nvars; k++) D[k][k] = absolute(D[k][k]);
      ierr = MatMult(nvars,DL,D,L);    CHECKERR(ierr);
      ierr = MatMult(nvars,modA,R,DL); CHECKERR(ierr);

      for (k=0; k<nvars; k++) {
        udiss[k] = 0; for (v=0; v<nvars; v++) udiss[k] += modA[k][v]*udiff[v];
      }
      
      for (v=0; v<nvars; v++) fI[nvars*p+v] = 0.5 * (fL[nvars*p+v]+fR[nvars*p+v]) - udiss[v];
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  for (k = 0; k < nvars; k++) {
    free(modA[k]);
    free(R[k]);
    free(D[k]);
    free(L[k]);
    free(DL[k]);
  }
  free(modA);
  free(R);
  free(D);
  free(L);
  free(DL);

  free(index_inter);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);

  return(0);
}

int Euler2DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      ierr    = 0,done,k;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *index_inter  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  int *bounds_inter = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  double **R, **D, **L;
  R     = (double**) calloc (nvars,sizeof(double*));
  D     = (double**) calloc (nvars,sizeof(double*));
  L     = (double**) calloc (nvars,sizeof(double*));
  for (k = 0; k < nvars; k++){
    R[k]    = (double*) calloc (nvars,sizeof(double));
    D[k]    = (double*) calloc (nvars,sizeof(double));
    L[k]    = (double*) calloc (nvars,sizeof(double));
  }

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double uavg[4], fcL[4], fcR[4], ucL[4], ucR[4], fc[4];

      /* Roe-Fixed upwinding scheme */

      ierr = Euler2DRoeAverage(&uavg[0],&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      ierr = Euler2DEigenvalues       (&uavg[0],D,param,dir); CHECKERR(ierr);
      ierr = Euler2DLeftEigenvectors  (&uavg[0],L,param,dir); CHECKERR(ierr);
      ierr = Euler2DRightEigenvectors (&uavg[0],R,param,dir); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      ierr = MatVecMult(nvars,&ucL[0],L,&uL[nvars*p]);
      ierr = MatVecMult(nvars,&ucR[0],L,&uR[nvars*p]);
      ierr = MatVecMult(nvars,&fcL[0],L,&fL[nvars*p]);
      ierr = MatVecMult(nvars,&fcR[0],L,&fR[nvars*p]);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        ierr = Euler2DEigenvalues(&uL[nvars*p],D,param,dir); CHECKERR(ierr);
        eigL = D[k][k];
        ierr = Euler2DEigenvalues(&uR[nvars*p],D,param,dir); CHECKERR(ierr);
        eigR = D[k][k];
        ierr = Euler2DEigenvalues(&uavg[0]    ,D,param,dir); CHECKERR(ierr);
        eigC = D[k][k];

        if ((eigL > 0) && (eigC > 0) && (eigR > 0)) {
          fc[k] = fcL[k];
        } else if ((eigL < 0) && (eigC < 0) && (eigR < 0)) {
          fc[k] = fcR[k];
        } else {
          double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      ierr = MatVecMult(nvars,&fI[nvars*p],R,&fc[0]); CHECKERR(ierr);
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  for (k = 0; k < nvars; k++) {
    free(R[k]);
    free(D[k]);
    free(L[k]);
  }
  free(R);
  free(D);
  free(L);

  free(index_inter);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);

  return(0);
}

int Euler2DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar           *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int             ierr    = 0,done,k;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *index_inter  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  int *bounds_inter = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  double **R, **D, **L;
  R     = (double**) calloc (nvars,sizeof(double*));
  D     = (double**) calloc (nvars,sizeof(double*));
  L     = (double**) calloc (nvars,sizeof(double*));
  for (k = 0; k < nvars; k++){
    R[k]    = (double*) calloc (nvars,sizeof(double));
    D[k]    = (double*) calloc (nvars,sizeof(double));
    L[k]    = (double*) calloc (nvars,sizeof(double));
  }

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double uavg[4], fcL[4], fcR[4], ucL[4], ucR[4], fc[4];

      /* Roe-Fixed upwinding scheme */

      ierr = Euler2DRoeAverage(&uavg[0],&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      ierr = Euler2DEigenvalues       (&uavg[0],D,param,dir); CHECKERR(ierr);
      ierr = Euler2DLeftEigenvectors  (&uavg[0],L,param,dir); CHECKERR(ierr);
      ierr = Euler2DRightEigenvectors (&uavg[0],R,param,dir); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      ierr = MatVecMult(nvars,&ucL[0],L,&uL[nvars*p]);
      ierr = MatVecMult(nvars,&ucR[0],L,&uR[nvars*p]);
      ierr = MatVecMult(nvars,&fcL[0],L,&fL[nvars*p]);
      ierr = MatVecMult(nvars,&fcR[0],L,&fR[nvars*p]);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        ierr = Euler2DEigenvalues(&uL[nvars*p],D,param,dir); CHECKERR(ierr);
        eigL = D[k][k];
        ierr = Euler2DEigenvalues(&uR[nvars*p],D,param,dir); CHECKERR(ierr);
        eigR = D[k][k];
        ierr = Euler2DEigenvalues(&uavg[0]    ,D,param,dir); CHECKERR(ierr);
        eigC = D[k][k];

        double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      ierr = MatVecMult(nvars,&fI[nvars*p],R,&fc[0]); CHECKERR(ierr);
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  for (k = 0; k < nvars; k++) {
    free(R[k]);
    free(D[k]);
    free(L[k]);
  }
  free(R);
  free(D);
  free(L);

  free(index_inter);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);

  return(0);
}
