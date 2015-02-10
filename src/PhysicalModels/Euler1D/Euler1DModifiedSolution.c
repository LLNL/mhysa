#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
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

int Euler1DModifiedSolution(double *uC,double *u,int d,void *s,void *m,double waqt)
{
  HyPar         *solver = (HyPar*)         s;
  Euler1D       *param  = (Euler1D*)       solver->physics;

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
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p); p *= _MODEL_NVARS_;
    _ArrayScaleCopy1D_((u+p),(1.0/param->grav_field[p/_MODEL_NVARS_]),(uC+p),_MODEL_NVARS_);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  return(0);
}
