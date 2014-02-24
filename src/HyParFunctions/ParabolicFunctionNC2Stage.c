#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ParabolicFunctionNC2Stage(double *par,double *u,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  double        *Func   = solver->fluxC;
  double        *Deriv1 = solver->Deriv1;
  double        *Deriv2 = solver->Deriv2;
  int           d, d1, d2, v, p, done;
  double        dxinv1, dxinv2;
  _DECLARE_IERR_;

  int     ndims   = solver->ndims;
  int     nvars   = solver->nvars;
  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  double  *dxinv  = solver->dxinv;

  if (!solver->HFunction) return(0); /* zero parabolic terms */

  int index[ndims];
  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  _ArraySetValue_(par,size*nvars,0.0);

  for (d1 = 0; d1 < ndims; d1++) {
    for (d2 = 0; d2 < ndims; d2++) {

      /* calculate the diffusion function */
      IERR solver->HFunction(Func,u,d1,d2,solver,t);                CHECKERR(ierr);
      IERR solver->FirstDerivativePar(Deriv1,Func  ,d1,solver,mpi); CHECKERR(ierr);
      IERR solver->FirstDerivativePar(Deriv2,Deriv1,d2,solver,mpi); CHECKERR(ierr);

      /* calculate the final term - second derivative of the diffusion function */
      done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        _ArrayIndex1D_(ndims,dim,index,ghosts,p);
        _GetCoordinate_(d1,index[d1],dim,ghosts,dxinv,dxinv1);
        _GetCoordinate_(d2,index[d2],dim,ghosts,dxinv,dxinv2);
        for (v=0; v<nvars; v++) par[nvars*p+v] = dxinv1*dxinv2 * Deriv2[nvars*p+v];
        _ArrayIncrementIndex_(ndims,dim,index,done);
      }

    }
  }

  return(0);
}
