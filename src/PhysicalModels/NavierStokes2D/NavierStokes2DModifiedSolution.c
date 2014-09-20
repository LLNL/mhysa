#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

/*

  Refer to 
  + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme 
    for the Gas Dynamics Equations Under Gravitational Fields", 
    Journal of Scientific Computing, 54, 2013, pp. 645-662,
    http://dx.doi.org/10.1007/s10915-012-9585-8

    Equation (3.8) on Page 651 on why this modification is needed.
*/

int NavierStokes2DModifiedSolution(double *uC,double *u,int d,void *s,void *m,double waqt)
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m; 
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  double inv_gamma_m1 = 1.0 / (param->gamma-1.0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, uvel, vvel, E, P;
    _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,E,P,param);
    uC[_MODEL_NVARS_*p+0] = u[_MODEL_NVARS_*p+0] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+1] = u[_MODEL_NVARS_*p+1] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+2] = u[_MODEL_NVARS_*p+2] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+3] = (P*inv_gamma_m1)*(1.0/param->grav_field_g[p]) + (0.5*rho*(uvel*uvel+vvel*vvel))*param->grav_field_f[p];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  return(0);
}
