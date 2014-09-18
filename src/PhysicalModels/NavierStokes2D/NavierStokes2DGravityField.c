#include <stdlib.h>
#include <math.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>
#include <mpivars.h>

int NavierStokes2DGravityField(void *s,void *m,double *u)
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;

  double  *f      = param->grav_field_f;
  double  *g      = param->grav_field_g;
  int     *dim    = solver->dim_local;
  int     ghosts  = solver->ghosts;
  int     index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_],
          offset[_MODEL_NDIMS_], d, done;

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_);
  for (d=0; d<_MODEL_NDIMS_; d++) bounds[d] += 2*ghosts;
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  double p0   = param->p0;
  double rho0 = param->rho0;
  double RT   = p0 / rho0;
  double gamma= param->gamma;
  double R    = param->R;
  double Cp   = gamma*R / (gamma-1.0);
  double T0   = p0 / (rho0 * R);
  double gx   = param->grav_x;
  double gy   = param->grav_y;
  
  /* set the value of the gravity field */
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  if (param->HB == 1) {
    while (!done) {
      int p;         _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
      double xcoord; _GetCoordinate_(_XDIR_,index[_XDIR_]-ghosts,dim,ghosts,solver->x,xcoord);
      double ycoord; _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->x,ycoord);
      f[p] = exp( (gx*xcoord+gy*ycoord)/RT);
      g[p] = exp(-(gx*xcoord+gy*ycoord)/RT);
      _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
    }
  } else if (param->HB == 2) {
    while (!done) {
      int p;         _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
      double xcoord; _GetCoordinate_(_XDIR_,index[_XDIR_]-ghosts,dim,ghosts,solver->x,xcoord);
      double ycoord; _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->x,ycoord);
      f[p] = raiseto((1.0-(gx*xcoord+gy*ycoord)/(Cp*T0)), (-1.0 /(gamma-1.0)));
      g[p] = raiseto((1.0-(gx*xcoord+gy*ycoord)/(Cp*T0)), (gamma/(gamma-1.0)));
      _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
    }
  }

  /* A sensible simulation will not specify peridic boundary conditions 
   * along a direction in which gravity acts.
   * Gravity will be zero along the dimension periodic BCs are specified,
   * so the value of the gravity potential will be the same along those
   * grid lines.
   * Thus, for both these cases, extrapolate the gravity field at the 
   * boundaries */
  int indexb[_MODEL_NDIMS_], indexi[_MODEL_NDIMS_];
  for (d = 0; d < _MODEL_NDIMS_; d++) {
    /* left boundary */
    if (!mpi->ip[d]) {
      _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
      _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = -ghosts;
      done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = ghosts-1-indexb[d];
        int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
        f[p1] = f[p2];
        g[p1] = g[p2];
        _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
      }
    }
    /* right boundary */
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
      _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = dim[d];
      done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = dim[d]-1-indexb[d];
        int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
        f[p1] = f[p2];
        g[p1] = g[p2];
        _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
      }
    }
  }

  return(0);
}
