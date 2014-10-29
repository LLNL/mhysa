#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ParabolicFunctionCons1Stage(double *par,double *u,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  double        *FluxI  = solver->fluxI; /* interface flux     array */
  double        *Func   = solver->fluxC; /* diffusion function array */
  int           d, v, i, done;
  _DECLARE_IERR_;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;

  int index[ndims], index1[ndims], index2[ndims], dim_interface[ndims];

  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  _ArraySetValue_(par,size*nvars,0.0);
  if (!solver->GFunction) return(0); /* zero parabolic term */

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[d]++;
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];

    /* evaluate cell-centered flux */
    IERR solver->GFunction(Func,u,d,solver,t); CHECKERR(ierr);
    /* compute interface fluxes */
    IERR solver->InterpolateInterfacesPar(FluxI,Func,d,solver,mpi); CHECKERR(ierr);

    /* calculate the second derivative */
    done = 0; _ArraySetValue_(index,ndims,0);
    int p, p1, p2;
    while (!done) {
      _ArrayCopy1D_(index,index1,ndims);
      _ArrayCopy1D_(index,index2,ndims); index2[d]++;
      _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p);
      _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
      _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
      for (v=0; v<nvars; v++) 
        par[nvars*p+v] +=  ((dxinv[offset+ghosts+index[d]] * dxinv[offset+ghosts+index[d]]) 
                          * (FluxI[nvars*p2+v] - FluxI[nvars*p1+v]));
      _ArrayIncrementIndex_(ndims,dim,index,done);
    }

    offset += dim[d] + 2*ghosts;
  }

  return(0);
}
