/*! @file BoundaryIntegral.c
    @author Debojyoti Ghosh
    @brief Compute the flux integral over the physical boundary.
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Computes the flux integral over the boundary. The local flux integral
    (on this processor) is computed for physical as well as MPI boundaries.
    The global boundary integral is computed by summing the local integrals
    over all the processors, since the contributions from the MPI boundaries
    cancel out.
*/
int BoundaryIntegral(
                      void *s, /*!< Solver object of type #HyPar */
                      void *m  /*!< MPI object of type #MPIVariables */
                    )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int d,v,k;

  double *local_integral  = (double*) calloc (nvars,sizeof(double));
  double *global_integral = (double*) calloc (nvars,sizeof(double));

  /* calculate the local boundary integral on each process */
  _ArraySetValue_(local_integral,nvars,0.0);
  for (d=0; d<ndims; d++) {
    for (v=0; v<nvars; v++) {
      double dxinv[ndims], dS = 1.0; 
      for (k=0; k<ndims; k++) { 
        /* uniform grid assumed */
        _GetCoordinate_(k,dim[k]/2,dim,ghosts,solver->dxinv,dxinv[k]); 
      }
      for (k=0; k<ndims; k++) if (k!=d) dS *= (1.0/dxinv[k]);
      local_integral[v] += (solver->StepBoundaryIntegral[(2*d+0)*nvars+v])*dS;
      local_integral[v] += (solver->StepBoundaryIntegral[(2*d+1)*nvars+v])*dS;
    }
  }

  /* add across process to calculate global boundary integral 
   * (internal (MPI) boundaries must cancel out)
   */
  IERR MPISum_double(global_integral,local_integral,nvars,&mpi->world); CHECKERR(ierr);

  /* add to the total boundary integral */
  _ArrayAXPY_(global_integral,1.0,solver->TotalBoundaryIntegral,nvars);

  free(local_integral);
  free(global_integral);
  return(0);
}
