#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

int ReconstructHyperbolic(double *fluxI,double *fluxC,double *u,int dir,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0, d;
  double        *fluxL  = NULL;
  double        *fluxR  = NULL;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  /* allocate arrays for left and right biased interface fluxes */
  int size = 1;
  for (d=0; d<ndims; d++) {
    if (d == dir) size  *= (dim[d]+1);
    else          size  *=  dim[d];
  }
  fluxL = (double*) calloc (size*nvars,sizeof(double));
  fluxR = (double*) calloc (size*nvars,sizeof(double));

  /* Interpolation -> to calculate left and right-biased interface flux */
  ierr = solver->InterpolateInterfacesHyp(fluxL,fluxC, 1,dir,solver,mpi); CHECKERR(ierr);
  ierr = solver->InterpolateInterfacesHyp(fluxR,fluxC,-1,dir,solver,mpi); CHECKERR(ierr);

  /* Upwind -> to calculate the final interface flux */
  ierr = solver->Upwind(fluxI,fluxL,fluxR,u,dir,solver,t); CHECKERR(ierr); 

  free(fluxL);
  free(fluxR);
  return(0);
}
