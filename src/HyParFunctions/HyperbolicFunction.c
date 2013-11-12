#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ReconstructHyperbolic (double*,double*,double*,int,void*,void*,double);

int HyperbolicFunction(double *hyp,double *u,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0, d, v, i, done;
  double        *FluxI  = solver->fluxI; /* interface flux     */
  double        *FluxC  = solver->fluxC; /* cell centered flux */

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;
  int     index[ndims], index1[ndims], index2[ndims], dim_interface[ndims];

  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  ierr = ArraySetValue_double(hyp,size*nvars,0.0); CHECKERR(ierr);
  if (!solver->FFunction) return(0); /* zero hyperbolic term */

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    ierr = ArrayCopy1D_int(dim,dim_interface,ndims); CHECKERR(ierr); dim_interface[d]++;
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];

    /* evaluate cell-centered flux */
    ierr = solver->FFunction(FluxC,u,d,solver,t); CHECKERR(ierr);
    /* compute interface fluxes */
    ierr = ReconstructHyperbolic(FluxI,FluxC,u,d,solver,mpi,t); CHECKERR(ierr);

    /* calculate the first derivative */
    done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
    while (!done) {
      ierr = ArrayCopy1D_int(index,index1,ndims); CHECKERR(ierr);
      ierr = ArrayCopy1D_int(index,index2,ndims); CHECKERR(ierr); index2[d]++;
      int p  = ArrayIndex1D(ndims,dim          ,index ,NULL,ghosts);
      int p1 = ArrayIndex1D(ndims,dim_interface,index1,NULL,0     );
      int p2 = ArrayIndex1D(ndims,dim_interface,index2,NULL,0     );
      for (v=0; v<nvars; v++) hyp[nvars*p+v] += dxinv[offset+ghosts+index[d]] 
                                              * (FluxI[nvars*p2+v]-FluxI[nvars*p1+v]);
      done = ArrayIncrementIndex(ndims,dim,index);
    }

    offset += dim[d] + 2*ghosts;
  }

  return(0);
}
