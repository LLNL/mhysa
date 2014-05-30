#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

int ReconstructHyperbolic(double *fluxI,double *fluxC,double *u,double *x,int dir,void *s,void *m,double t,int LimFlag)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           d;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  /* allocate arrays for left and right biased interface fluxes */
  int size = 1;
  for (d=0; d<ndims; d++) {
    if (d == dir) size  *= (dim[d]+1);
    else          size  *=  dim[d];
  }
  size *= nvars;
  double *uL     = solver->uL;
  double *uR     = solver->uR;
  double *fluxL  = solver->fL;
  double *fluxR  = solver->fR;

  /* 
    precalculate the non-linear interpolation coefficients if required 
    else reuse the weights previously calculated
  */
  if (LimFlag) IERR solver->SetInterpLimiterVar(fluxC,u,x,dir,solver,mpi);

  /* Interpolation -> to calculate left and right-biased interface flux and state variable*/
  IERR solver->InterpolateInterfacesHyp(uL   ,u    ,u,x, 1,dir,solver,mpi); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(uR   ,u    ,u,x,-1,dir,solver,mpi); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(fluxL,fluxC,u,x, 1,dir,solver,mpi); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(fluxR,fluxC,u,x,-1,dir,solver,mpi); CHECKERR(ierr);

  /* Upwind -> to calculate the final interface flux */
  IERR solver->Upwind(fluxI,fluxL,fluxR,uL,uR,u,dir,solver,t); CHECKERR(ierr); 

  return(0);
}
