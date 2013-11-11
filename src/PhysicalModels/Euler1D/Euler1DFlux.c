#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

inline int Euler1DGetFlowVar (double*,double*,double*,double*,double*,void*);
inline int Euler1DSetFlux    (double*,double ,double ,double ,double ,void*);

int Euler1DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;
  int       ierr    = 0, i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  ierr = ArrayCopy1D_int(dim,bounds,ndims); CHECKERR(ierr);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  ierr = ArraySetValue_int(offset,ndims,-ghosts); CHECKERR(ierr);

  int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,offset,ghosts);
    double rho, v, e, P;
    ierr = Euler1DGetFlowVar(&u[nvars*p],&rho,&v,&e,&P,param); CHECKERR(ierr);
    ierr = Euler1DSetFlux   (&f[nvars*p],rho ,v ,e ,P ,param); CHECKERR(ierr);
    done = ArrayIncrementIndex(ndims,bounds,index);
  }

  return(0);
}
