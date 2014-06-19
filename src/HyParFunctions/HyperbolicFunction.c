#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ReconstructHyperbolic (double*,double*,double*,double*,int,void*,void*,double,int);

int HyperbolicFunction(double *hyp,double *u,void *s,void *m,double t,int LimFlag,
                       int(*FluxFunction)(double*,double*,int,void*,double))
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

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[d]++;
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];

    /* evaluate cell-centered flux */
    IERR FluxFunction(FluxC,u,d,solver,t); CHECKERR(ierr);
    /* compute interface fluxes */
    IERR ReconstructHyperbolic(FluxI,FluxC,u,x+offset,d,solver,mpi,t,LimFlag); CHECKERR(ierr);

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
