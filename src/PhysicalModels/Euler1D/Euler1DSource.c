#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

/* function to calculate the gravitational source term
 * Reference: Xing, Shu, "High Order Well-Balanced WENO
 * Scheme for the Gas Dynamics Equations Under Gravitational
 * Fields", J. Sci. Comput., 54, 2013, pp. 645--662
 * http://dx.doi.org/10.1007/s10915-012-9585-8
*/

static int Euler1DSourceFunction(double*,double*,double*,void*,double);

int Euler1DSource(double *source,double *u,void *s,void *m,double t)
{
  HyPar         *solver = (HyPar* ) s;
  MPIVariables  *mpi = (MPIVariables*) m;
  Euler1D       *param  = (Euler1D*) solver->physics;

  if (param->grav == 0.0)  return(0); /* no gravitational forces */

  int     i, d, v, done, p, p1, p2;
  double  *SourceI = solver->fluxI; /* interace source term       */
  double  *SourceC = solver->fluxC; /* cell-centered source term  */
  double  *SourceL = solver->fL;
  double  *SourceR = solver->fR;

  int     ndims   = solver->ndims;
  int     nvars   = solver->nvars;
  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  double  *x      = solver->x;
  double  *dxinv  = solver->dxinv;
  int     index[ndims],index1[ndims],index2[ndims],dim_interface[ndims];

  /* set interface dimensions */
  _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_XDIR_]++;
  /* calculate the split source function exp(-phi/RT) */
  IERR Euler1DSourceFunction(SourceC,u,x,solver,t); CHECKERR(ierr);
  /* calculate the left and right interface source terms */
  IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_XDIR_,solver,mpi); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_XDIR_,solver,mpi); CHECKERR(ierr);
  /* calculate the final interface source term */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayIndex1D_(ndims,dim_interface,index,0,p);
    for (v=0; v<_MODEL_NVARS_;v++) 
      SourceI[_MODEL_NVARS_*p+v] = 0.5 * (SourceL[_MODEL_NVARS_*p+v] + SourceR[_MODEL_NVARS_*p+v]);
  }
  /* calculate the final cell-centered source term */
  while (!done) {
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_XDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index,0     ,p2);
    double dx_inverse;   _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dx_inverse);
    double rho, vel, e, P; _Euler1DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vel,e,P,param);
    double term[_MODEL_NVARS_] = {0.0, P, P*vel};
    for (v=0; v<_MODEL_NVARS_;v++) {
      source[_MODEL_NVARS_*p+v] += (  (term[v]*(1.0/SourceC[_MODEL_NVARS_*p+v])) 
                                    * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse );
    }
    
  }

  return(0);
}

int Euler1DSourceFunction(double *f,double *u,double *x,void *s,double t)
{
  HyPar   *solver = (HyPar* ) s;
  Euler1D *param  = (Euler1D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i;
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, v, e, P, xcoord;
    _GetCoordinate_(_XDIR_,index[_XDIR_]-ghosts,dim,ghosts,x,xcoord);
    _Euler1DGetFlowVar_((u+_MODEL_NVARS_*p),rho,v,e,P,param);
    _Euler1DSetSource_((f+_MODEL_NVARS_*p),rho,P,param,xcoord);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
