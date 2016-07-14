/*! @file InitialSolution.c
    @author Debojyoti Ghosh
    @brief Read in initial solution from file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <io.h>
#include <mpivars.h>
#include <hypar.h>

int VolumeIntegral(double*,double*,void*,void*);

/*! Read in initial solution from file, and compute grid spacing 
    and volume integral of the initial solution */
int InitialSolution(void *s, void *m)
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           flag, d, i, offset, ghosts = solver->ghosts;
  _DECLARE_IERR_;

  IERR ReadArray(solver->ndims,solver->nvars,solver->dim_global,solver->dim_local,
                 solver->ghosts,solver,mpi,solver->x,solver->u,"initial",&flag);
  if (!flag) {
    fprintf(stderr,"Error: initial solution file not found.\n");
    return(1);
  }
  CHECKERR(ierr);

  /* exchange MPI-boundary values of u between processors */
  MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,solver->ghosts,mpi,solver->u);

  /* calculate dxinv */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    for (i = 0; i < solver->dim_local[d]; i++) 
      solver->dxinv[i+offset+ghosts] = 2.0/(solver->x[i+1+offset+ghosts]-solver->x[i-1+offset+ghosts]);
    offset += (solver->dim_local[d] + 2*ghosts);
  }

  /* exchange MPI-boundary values of dxinv between processors */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->dxinv[offset],solver->dim_local[d],
                                   ghosts,d,solver->ndims); CHECKERR(ierr);
    offset  += (solver->dim_local[d] + 2*ghosts);
  }

  /* fill in ghost values of dxinv at physical boundaries by extrapolation */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *dxinv = &solver->dxinv[offset];
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) dxinv[i] = dxinv[ghosts];
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) dxinv[i] = dxinv[dim[d]+ghosts-1];
    }
    offset  += (dim[d] + 2*ghosts);
  }

  /* calculate volume integral of the initial solution */
  IERR VolumeIntegral(solver->VolumeIntegralInitial,solver->u,solver,mpi); CHECKERR(ierr);
  if (!mpi->rank) {
    printf("Volume integral of the initial solution:\n");
    for (d=0; d<solver->nvars; d++) printf("%2d:  %1.16E\n",d,solver->VolumeIntegralInitial[d]);
  }
  /* Set initial total boundary flux integral to zero */
  _ArraySetValue_(solver->TotalBoundaryIntegral,solver->nvars,0);

  return(0); 
}
