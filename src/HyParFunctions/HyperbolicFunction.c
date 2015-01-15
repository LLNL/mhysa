#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ReconstructHyperbolic (double*,double*,double*,double*,int,void*,void*,double,int,
                                  int(*)(double*,double*,double*,double*,double*,double*,int,void*,double));
static int DefaultUpwinding      (double*,double*,double*,double*,double*,double*,int,void*,double);

int HyperbolicFunction(double *hyp,double *u,void *s,void *m,double t,int LimFlag,
                       int(*FluxFunction)(double*,double*,int,void*,double),
                       int(*UpwindFunction)(double*,double*,double*,double*,double*,double*,int,void*,double))
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           d, v, i, done;
  double        *FluxI  = solver->fluxI; /* interface flux     */
  double        *FluxC  = solver->fluxC; /* cell centered flux */
  _DECLARE_IERR_;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *x     = solver->x;
  double  *dxinv = solver->dxinv;
  int     index[ndims], index1[ndims], index2[ndims], dim_interface[ndims];

  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  LimFlag = (LimFlag && solver->flag_nonlinearinterp && solver->SetInterpLimiterVar);

  _ArraySetValue_(hyp,size*nvars,0.0);
  _ArraySetValue_(solver->StageBoundaryIntegral,2*ndims*nvars,0.0);
  if (!FluxFunction) return(0); /* zero hyperbolic term */
  solver->count_hyp++;

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[d]++;
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];

    /* evaluate cell-centered flux */
    IERR FluxFunction(FluxC,u,d,solver,t); CHECKERR(ierr);
    /* compute interface fluxes */
    IERR ReconstructHyperbolic(FluxI,FluxC,u,x+offset,d,solver,mpi,t,LimFlag,UpwindFunction); 
    CHECKERR(ierr);

    /* calculate the first derivative */
    done = 0; _ArraySetValue_(index,ndims,0);
    int p, p1, p2;
    while (!done) {
      _ArrayCopy1D_(index,index1,ndims);
      _ArrayCopy1D_(index,index2,ndims); index2[d]++;
      _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p);
      _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
      _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
      for (v=0; v<nvars; v++) hyp[nvars*p+v] += dxinv[offset+ghosts+index[d]] 
                                              * (FluxI[nvars*p2+v]-FluxI[nvars*p1+v]);
      /* boundary flux integral */
      if (index[d] == 0) 
        for (v=0; v<nvars; v++) solver->StageBoundaryIntegral[(2*d+0)*nvars+v] -= FluxI[nvars*p1+v];
      if (index[d] == dim[d]-1) 
        for (v=0; v<nvars; v++) solver->StageBoundaryIntegral[(2*d+1)*nvars+v] += FluxI[nvars*p2+v];

      _ArrayIncrementIndex_(ndims,dim,index,done);
    }

    offset += dim[d] + 2*ghosts;
  }

  return(0);
}

int ReconstructHyperbolic(double *fluxI,double *fluxC,double *u,double *x,int dir,void *s,void *m,double t,int LimFlag,
                          int(*UpwindFunction)(double*,double*,double*,double*,double*,double*,int,void*,double))
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           d;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  double *uC     = NULL;
  double *uL     = solver->uL;
  double *uR     = solver->uR;
  double *fluxL  = solver->fL;
  double *fluxR  = solver->fR;

  /* 
    precalculate the non-linear interpolation coefficients if required 
    else reuse the weights previously calculated
  */
  if (LimFlag) IERR solver->SetInterpLimiterVar(fluxC,u,x,dir,solver,mpi);

  /* if defined, calculate the modified u-function to be used for upwinding
     e.g.: used in well-balanced schemes for Euler/Navier-Stokes with gravity
     otherwise, just copy u to uC */
  if (solver->UFunction) {
    uC = solver->uC;
    IERR solver->UFunction(uC,u,d,solver,mpi,t); CHECKERR(ierr);
  } else uC = u;

  /* Interpolation -> to calculate left and right-biased interface flux and state variable*/
  IERR solver->InterpolateInterfacesHyp(uL   ,uC   ,u,x, 1,dir,solver,mpi,1); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(uR   ,uC   ,u,x,-1,dir,solver,mpi,1); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(fluxL,fluxC,u,x, 1,dir,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(fluxR,fluxC,u,x,-1,dir,solver,mpi,0); CHECKERR(ierr);

  /* Upwind -> to calculate the final interface flux */
  if (UpwindFunction) { IERR UpwindFunction   (fluxI,fluxL,fluxR,uL  ,uR  ,u   ,dir,solver,t); CHECKERR(ierr); }
  else                { IERR DefaultUpwinding (fluxI,fluxL,fluxR,NULL,NULL,NULL,dir,solver,t); CHECKERR(ierr); }

  return(0);
}

int DefaultUpwinding(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar *solver = (HyPar*)    s;
  int   done;

  int *dim  = solver->dim_local;
  int ndims = solver->ndims;
  int nvars = solver->nvars;

  int bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; int index_outer[ndims], index_inter[ndims];
  _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int v; for (v=0; v<nvars; v++) fI[nvars*p+v] = 0.5 * (fL[nvars*p+v]+fR[nvars*p+v]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
