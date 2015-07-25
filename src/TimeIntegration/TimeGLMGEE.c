#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeGLMGEE(void *ts)
{
  TimeIntegration   *TS     = (TimeIntegration*) ts;
  HyPar             *solver = (HyPar*)           TS->solver;
  MPIVariables      *mpi    = (MPIVariables*)    TS->mpi;
  GLMGEEParameters  *params = (GLMGEEParameters*)solver->msti;
  int               i, j;
  _DECLARE_IERR_;

  int    s      = params->nstages;
  int    r      = params->r;
  int    ndims  = solver->ndims;
  int    nvars  = solver->nvars;
  int    size   = nvars * solver->npoints_local_wghosts;
  double dt     = TS->dt;
  double *A=params->A, *B=params->B, *C=params->C, *D=params->D, *c=params->c,
         **U = TS->U, **Udot = TS->Udot, **Uaux = &TS->U[r];

  /* Calculate stage values */
  for (j=0; j<s; j++) {

    double stagetime = TS->waqt + c[j]*dt;

    _ArrayScaleCopy1D_(solver->u,C[j*r+0],U[0],size);
    for (i=1;i<r;i++) _ArrayAXPY_(Uaux[i-1],C[j*r+i]   ,U[0],size);
    for (i=0;i<j;i++) _ArrayAXPY_(Udot[i]  ,dt*A[j*s+i],U[0],size); 

    if (solver->PreStage) { IERR solver->PreStage(j,U,solver,mpi,stagetime); CHECKERR(ierr); }
    IERR TS->RHSFunction(Udot[j],U[0],solver,mpi,stagetime);
    if (solver->PostStage) { IERR solver->PostStage(U[j],solver,mpi,stagetime); CHECKERR(ierr); }

    _ArraySetValue_(TS->BoundaryFlux[j],2*ndims*nvars,0.0);
    _ArrayCopy1D_(solver->StageBoundaryIntegral,TS->BoundaryFlux[j],2*ndims*nvars);
  }

  /* Step completion */
  for (j=0; j<r; j++) {
    _ArrayScaleCopy1D_(solver->u,D[j*r+0],U[j],size);
    for (i=1; i<r; i++) _ArrayAXPY_(Uaux[i-1],D[j*r+i]   ,U[j],size);
    for (i=0; i<s; i++) _ArrayAXPY_(Udot[i]  ,dt*B[j*s+i],U[j],size);
  }

  for (i=0; i<s; i++) _ArrayAXPY_(TS->BoundaryFlux[i],dt*B[0*s+i],solver->StepBoundaryIntegral,
                                  2*ndims*nvars);

  _ArrayCopy1D_(U[0],solver->u,size);
  for (i=1; i<r; i++) _ArrayCopy1D_(U[i],Uaux[i-1],size);

  return(0);
}

