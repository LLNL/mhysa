/*! @file VolumeIntegral.c
    @author Debojyoti Ghosh
    @brief Compute the volume integral of the solution
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the volume integral of the solution.
*/
int VolumeIntegral(
                    double  *VolumeIntegral,  /*!< The computed volume integral */
                    double  *u,               /*!< Solution */
                    void    *s,               /*!< Solver object of type #HyPar */
                    void    *m                /*!< MPI object of type #MPIVariables */
                  )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int d,v;

  /* calculate local volume integral of the solution */
  int index[ndims], done = 0;
  double *local_integral = (double*) calloc (nvars,sizeof(double));
  _ArraySetValue_(local_integral,nvars,0.0);
  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dxinv[ndims];
    for (d=0; d<ndims; d++) { _GetCoordinate_(d,index[d],dim,ghosts,solver->dxinv,dxinv[d]); }
    double dV = 1.0; for (d=0; d<ndims; d++) dV *= (1.0/dxinv[d]);
    for (v=0; v<nvars; v++) local_integral[v] += (u[p*nvars+v]*dV);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  /* sum over all processors to get global integral of the solution */
  IERR MPISum_double(VolumeIntegral,local_integral,nvars,&mpi->world); CHECKERR(ierr);
  free(local_integral);

  return(0);
}
