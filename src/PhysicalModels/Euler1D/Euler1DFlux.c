#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DGetFlowVar (double*,double*,double*,double*,double*,void*);
int Euler1DSetFlux    (double*,double ,double ,double ,double ,void*);

int Euler1DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;
  int       i;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, v, e, P;
    IERR Euler1DGetFlowVar(&u[nvars*p],&rho,&v,&e,&P,param); CHECKERR(ierr);
    IERR Euler1DSetFlux   (&f[nvars*p],rho ,v ,e ,P ,param); CHECKERR(ierr);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
