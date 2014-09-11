#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

/* function to calculate the gravitational source term
 * Reference: Xing, Shu, "High Order Well-Balanced WENO
 * Scheme for the Gas Dynamics Equations Under Gravitational
 * Fields", J. Sci. Comput., 54, 2013, pp. 645--662
 * http://dx.doi.org/10.1007/s10915-012-9585-8
*/

static int NavierStokes2DSourceFunction (double*,double*,double*,void*,void*,double,int);
static int NavierStokes2DSourceUpwind   (double*,double*,double*,double*,int,void*,double);

int NavierStokes2DSource(double *source,double *u,void *s,void *m,double t)
{
  HyPar           *solver = (HyPar* )         s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;

  if ((param->grav_x == 0.0) && (param->grav_y == 0.0))  
    return(0); /* no gravitational forces */

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
  double  RT      =  param->p0 / param->rho0;
  int     index[ndims],index1[ndims],index2[ndims],dim_interface[ndims];

  /* Along X-direction */
  if (param->grav_x != 0.0) {
    /* set interface dimensions */
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_XDIR_]++;
    /* calculate the split source function exp(-phi/RT) */
    IERR NavierStokes2DSourceFunction(SourceC,u,x,solver,mpi,t,_XDIR_); CHECKERR(ierr);
    /* calculate the left and right interface source terms */
    IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
    IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
    /* calculate the final interface source term */
    IERR NavierStokes2DSourceUpwind(SourceI,SourceL,SourceR,u,_XDIR_,solver,t);
    /* calculate the final cell-centered source term */
    done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      _ArrayCopy1D_(index,index1,ndims);
      _ArrayCopy1D_(index,index2,ndims); index2[_XDIR_]++;
      _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
      _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
      _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);

      double dx_inverse; _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dx_inverse);
      double rho, uvel, vvel, e, P; _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,e,P,param);
      double term[_MODEL_NVARS_] = {0.0, rho*RT, 0.0, rho*RT*uvel};
      for (v=0; v<_MODEL_NVARS_; v++) {
        source[_MODEL_NVARS_*p+v] += (  (term[v]*(1.0/param->grav_field[p])) 
                                      * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse );
      }
      _ArrayIncrementIndex_(ndims,dim,index,done);
    }
  }

  /* Along Y-direction */
  if (param->grav_y != 0.0) {
    /* set interface dimensions */
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_YDIR_]++;
    /* calculate the split source function exp(-phi/RT) */
    IERR NavierStokes2DSourceFunction(SourceC,u,x,solver,mpi,t,_YDIR_); CHECKERR(ierr);
    /* calculate the left and right interface source terms */
    IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_YDIR_,solver,mpi,0); CHECKERR(ierr);
    IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_YDIR_,solver,mpi,0); CHECKERR(ierr);
    /* calculate the final interface source term */
    IERR NavierStokes2DSourceUpwind(SourceI,SourceL,SourceR,u,_YDIR_,solver,t);
    /* calculate the final cell-centered source term */
    done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      _ArrayCopy1D_(index,index1,ndims);
      _ArrayCopy1D_(index,index2,ndims); index2[_YDIR_]++;
      _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
      _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
      _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);

      double dy_inverse; _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,dxinv,dy_inverse);
      double rho, uvel, vvel, e, P; _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,e,P,param);
      double term[_MODEL_NVARS_] = {0.0, 0.0, rho*RT, rho*RT*vvel};
      for (v=0; v<_MODEL_NVARS_; v++) {
        source[_MODEL_NVARS_*p+v] += (  (term[v]*(1.0/param->grav_field[p])) 
                                      * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dy_inverse );
      }
      _ArrayIncrementIndex_(ndims,dim,index,done);
    }
  }

  return(0);
}

int NavierStokes2DSourceFunction(double *f,double *u,double *x,void *s,void *m,double t,int dir)
{
  HyPar           *solver = (HyPar* )         s;
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
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    (f+_MODEL_NVARS_*p)[0] = 0.0;
    (f+_MODEL_NVARS_*p)[1] = param->grav_field[p] * (dir == _XDIR_);
    (f+_MODEL_NVARS_*p)[2] = param->grav_field[p] * (dir == _YDIR_);
    (f+_MODEL_NVARS_*p)[3] = param->grav_field[p];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

int NavierStokes2DSourceUpwind(double *fI,double *fL,double *fR,double *u,int dir,void *s,double t)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (k = 0; k < _MODEL_NVARS_; k++) 
        (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
