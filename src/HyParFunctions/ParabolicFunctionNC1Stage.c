#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ParabolicFunctionNC1Stage(double *par,double *u,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0, d, v, i, done;
  double        *Func   = NULL;
  double        *Deriv2 = NULL;

  int     ndims   = solver->ndims;
  int     nvars   = solver->nvars;
  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  double  *dxinv  = solver->dxinv;

  int *index = (int*) calloc (ndims,sizeof(int));

  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  ierr = ArraySetValue_double(par,size*nvars,0.0); CHECKERR(ierr);
  if (!solver->GFunction) return(0); /* zero parabolic terms */

  int offset = 0;
  for (d = 0; d < ndims; d++) {

    /* allocate array for the diffusion function */
    Func = (double*) calloc (size*nvars,sizeof(double));
    /* calculate the diffusion function */
    ierr = solver->GFunction(Func,u,d,solver,t); CHECKERR(ierr);

    /* allocate array for the second differences and calculate it */
    int size_deriv = 1; for (i=0; i<ndims; i++) size_deriv *= dim[i];
    Deriv2  = (double*) calloc (size_deriv*nvars,sizeof(double));
    ierr    = solver->SecondDerivativePar(Deriv2,Func,d,solver,mpi); CHECKERR(ierr);

    /* calculate the final term - second derivative of the diffusion function */
    done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
    while (!done) {
      int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
      int q = ArrayIndex1D(ndims,dim,index,NULL,0     );
      for (v=0; v<nvars; v++)
        par[nvars*p+v] = dxinv[offset+index[d]]*dxinv[offset+index[d]] * Deriv2[nvars*q+v];
      done = ArrayIncrementIndex(ndims,dim,index);
    }

    /* free allocated arrays */
    free(Func);
    free(Deriv2);

    offset += dim[d];
  }

  free(index);
  return(0);
}
