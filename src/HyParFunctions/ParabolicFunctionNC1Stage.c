#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ParabolicFunctionNC1Stage(double *par,double *u,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  double        *Func   = solver->fluxC;
  double        *Deriv2 = solver->fluxI;
  int           d, v, i, done;
  _DECLARE_IERR_;

  int     ndims   = solver->ndims;
  int     nvars   = solver->nvars;
  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  double  *dxinv  = solver->dxinv;

  if (!solver->GFunction) return(0); /* zero parabolic terms */

  int index[ndims];
  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  _ArraySetValue_(par,size*nvars,0.0);

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    int size_deriv = 1; for (i=0; i<ndims; i++) size_deriv *= dim[i];

    /* calculate the diffusion function */
    IERR solver->GFunction(Func,u,d,solver,t); CHECKERR(ierr);
    IERR solver->SecondDerivativePar(Deriv2,Func,d,solver,mpi); CHECKERR(ierr);

    /* calculate the final term - second derivative of the diffusion function */
    done = 0; _ArraySetValue_(index,ndims,0);
    int p, q;
    while (!done) {
      _ArrayIndex1D_(ndims,dim,index,ghosts,p);
      _ArrayIndex1D_(ndims,dim,index,0     ,q);
      for (v=0; v<nvars; v++)
        par[nvars*p+v] =    dxinv[offset+ghosts+index[d]]*dxinv[offset+ghosts+index[d]] 
                          * Deriv2[nvars*q+v];
      _ArrayIncrementIndex_(ndims,dim,index,done);
    }

    offset += dim[d] + 2*ghosts;
  }

  return(0);
}
