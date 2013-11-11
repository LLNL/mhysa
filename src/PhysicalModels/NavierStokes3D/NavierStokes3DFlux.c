#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

inline int NavierStokes3DGetFlowVar (double*,double*,double*,double*,double*,double*,double*,void*);
inline int NavierStokes3DSetFlux    (double*,double ,double ,double ,double ,double ,double ,void*,int);

int NavierStokes3DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;
  int               ierr    = 0, i;

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
    double rho, vx, vy, vz, e, P;
    ierr = NavierStokes3DGetFlowVar(&u[nvars*p],&rho,&vx,&vy,&vz,&e,&P,param);     CHECKERR(ierr);
    ierr = NavierStokes3DSetFlux   (&f[nvars*p],rho ,vx ,vy ,vz ,e ,P ,param,dir); CHECKERR(ierr);
    done = ArrayIncrementIndex(ndims,bounds,index);
  }
  
  return(0);
}
