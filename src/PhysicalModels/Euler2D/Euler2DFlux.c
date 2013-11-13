#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

int Euler2DGetFlowVar (double*,double*,double*,double*,double*,double*,void*);
int Euler2DSetFlux    (double*,double ,double ,double ,double ,double ,void*,int);

int Euler2DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar   *solver = (HyPar*)   s;
  Euler2D *param  = (Euler2D*) solver->physics;
  int     i;
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
    double rho, vx, vy, e, P;
    IERR Euler2DGetFlowVar(&u[nvars*p],&rho,&vx,&vy,&e,&P,param);     CHECKERR(ierr);
    IERR Euler2DSetFlux   (&f[nvars*p],rho ,vx ,vy ,e ,P ,param,dir); CHECKERR(ierr);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
