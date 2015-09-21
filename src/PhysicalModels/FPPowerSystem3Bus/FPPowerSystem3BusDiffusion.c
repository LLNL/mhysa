/*! @file FPPowerSystem3BusDiffusion.c
    @author Debojyoti Ghosh
    @brief Compute the dissipative term for the #FPPowerSystem3Bus system.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

int FPPowerSystem3BusDissipationFunction(int,int,void*,double,double*);

/*! Compute the dissipation term for the #FPPowerSystem3Bus system */
int FPPowerSystem3BusDiffusion(
                                double  *f,   /*!< Array to hold the computed dissipation term vector (same layout as u) */
                                double  *u,   /*!< Array with the solution vector */
                                int     dir1, /*!< First spatial dimension for the dissipation term being computed */
                                int     dir2, /*!< Second spatial dimension for the dissipation term being computed */
                                void    *s,   /*!< Solver object of type #HyPar */
                                double  t     /*!< Current simulation time */
                              )
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->physics;
  int               i, v;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  static int index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];
  static double dissipation[_MODEL_NDIMS_*_MODEL_NDIMS_]; 

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_);
  for (i=0; i<_MODEL_NDIMS_; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    FPPowerSystem3BusDissipationFunction(dir1,dir2,params,t,dissipation);
    for (v = 0; v < _MODEL_NVARS_; v++) f[_MODEL_NVARS_*p+v] = dissipation[dir1*_MODEL_NDIMS_+dir2] * u[_MODEL_NVARS_*p+v];
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}
