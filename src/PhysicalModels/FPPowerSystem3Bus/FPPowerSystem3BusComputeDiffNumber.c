/*! @file FPPowerSystem3BusComputeDiffNumber.c
    @author Debojyoti Ghosh
    @brief Function to compute the maximum diffusion number.
*/

#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <mpivars.h>
#include <hypar.h>

int FPPowerSystem3BusDissipationFunction(int,int,void*,double,double*);

/*! Computes the maximum diffusion number over the domain. Note that the diffusion
    is computed over the local domain on this processor only.
*/
double FPPowerSystem3BusComputeDiffNumber(
                                          void    *s, /*!< Solver object of type #HyPar */
                                          void    *m, /*!< MPI object of type #MPIVariables */
                                          double  dt, /*!< Time step size for which to compute the CFL */
                                          double  t   /*!< Time */
                                         )
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->physics;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_diff = 0;
  int     index[ndims];
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    double dxinv[ndims],dissp[ndims*ndims];
    _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv[0]);
    _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dxinv[1]);
    _GetCoordinate_(2,index[2],dim,ghosts,solver->dxinv,dxinv[2]);
    _GetCoordinate_(3,index[3],dim,ghosts,solver->dxinv,dxinv[3]);
    FPPowerSystem3BusDissipationFunction(0,0,params,t,dissp);

    double local_diff[ndims];
    local_diff[0] = absolute(dissp[0*ndims+0]) * dt * dxinv[0] * dxinv[0];
    local_diff[1] = absolute(dissp[1*ndims+1]) * dt * dxinv[1] * dxinv[1];
    local_diff[2] = absolute(dissp[2*ndims+2]) * dt * dxinv[2] * dxinv[2];
    local_diff[3] = absolute(dissp[3*ndims+3]) * dt * dxinv[3] * dxinv[3];

    if (local_diff[0] > max_diff) max_diff = local_diff[0];
    if (local_diff[1] > max_diff) max_diff = local_diff[1];
    if (local_diff[2] > max_diff) max_diff = local_diff[2];
    if (local_diff[3] > max_diff) max_diff = local_diff[3];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_diff);
}

