/*! @file Initialize.c
    @author Debojyoti Ghosh
    @brief Initialization function
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Initialization function called at the beginning of a simulation. This function
    does the following:
    + allocates memory for MPI related arrays
    + initializes the values for MPI variables
    + creates sub-communicators and communication groups
    + allocates memory for arrays to store solution, right-hand-side, 
      flux, and other working vectors.
    + initializes function counters to zero
*/
int Initialize(
                void *s, /*!< Solver object of type #HyPar */
                void *m  /*!< MPI object of type #MPIVariables */
              )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i,d;

  /* allocations */
  mpi->ip           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->is           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->ie           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->bcperiodic   = (int*) calloc (solver->ndims,sizeof(int));
  solver->dim_local = (int*) calloc (solver->ndims,sizeof(int));
  solver->isPeriodic= (int*) calloc (solver->ndims,sizeof(int));

#ifndef serial
  _DECLARE_IERR_;

  /* Domain partitioning */
  if (!mpi->rank) printf("Partitioning domain.\n");
  int total_proc = 1;
  for (i=0; i<solver->ndims; i++) total_proc *= mpi->iproc[i];
  if (mpi->nproc != total_proc) {
    fprintf(stderr,"Error on rank %d: total number of processes is not consistent ", mpi->rank);
    fprintf(stderr,"with number of processes along each dimension.\n");
    fprintf(stderr,"mpiexec was called with %d processes, ",mpi->nproc);
    fprintf(stderr,"total number of processes from \"solver.inp\" is %d.\n", total_proc);
    return(1);
  }

  /* calculate ndims-D rank of each process (ip[]) from rank in MPI_COMM_WORLD */
  IERR MPIRanknD(solver->ndims,mpi->rank,mpi->iproc,mpi->ip); CHECKERR(ierr);

  /* calculate local domain sizes along each dimension */
  for (i=0; i<solver->ndims; i++) 
    solver->dim_local[i] = MPIPartition1D(solver->dim_global[i],mpi->iproc[i],mpi->ip[i]);

  /* calculate local domain limits in terms of global domain */
  IERR MPILocalDomainLimits(solver->ndims,mpi->rank,mpi,solver->dim_global,mpi->is,mpi->ie);
  CHECKERR(ierr);

  /* create sub-communicators for parallel computations along grid lines in each dimension */
  IERR MPICreateCommunicators(solver->ndims,mpi); CHECKERR(ierr);

  /* initialize periodic BC flags to zero */
  for (i=0; i<solver->ndims; i++) mpi->bcperiodic[i] = 0;

  /* create communication groups */
  IERR MPICreateIOGroups(mpi); CHECKERR(ierr);

#else

  for (i=0; i<solver->ndims; i++) {
    mpi->ip[i]            = 0;
    solver->dim_local[i]  = solver->dim_global[i];
    mpi->iproc[i]         = 1;
    mpi->is[i]            = 0;
    mpi->ie[i]            = solver->dim_local[i];
    mpi->bcperiodic[i]    = 0;
  }

#endif

  solver->npoints_global = solver->npoints_local = solver->npoints_local_wghosts = 1;
  for (i=0; i<solver->ndims; i++) solver->npoints_global *= solver->dim_global[i];
  for (i=0; i<solver->ndims; i++) solver->npoints_local  *= solver->dim_local [i];
  for (i=0; i<solver->ndims; i++) 
    solver->npoints_local_wghosts *= (solver->dim_local[i]+2*solver->ghosts);

  /* Allocations */
  if (!mpi->rank) printf("Allocating data arrays.\n");
  solver->index = (int*) calloc (solver->ndims,sizeof(int));
  solver->stride_with_ghosts    = (int*) calloc (solver->ndims,sizeof(int));
  solver->stride_without_ghosts = (int*) calloc (solver->ndims,sizeof(int));
  int accu1 = 1, accu2 = 1;
  for (i=0; i<solver->ndims; i++) {
    solver->stride_with_ghosts[i]    = accu1;
    solver->stride_without_ghosts[i] = accu2;
    accu1 *= (solver->dim_local[i]+2*solver->ghosts);
    accu2 *=  solver->dim_local[i];
  }
  int size;
  /* state variable */
  size = 1;
  for (i=0; i<solver->ndims; i++) size *= (solver->dim_local[i]+2*solver->ghosts);
  solver->u       = (double*) calloc (solver->nvars*size,sizeof(double));
#ifdef with_petsc
  if (solver->use_petscTS) {
    solver->u0      = (double*) calloc (solver->nvars*size,sizeof(double));
    solver->uref    = (double*) calloc (solver->nvars*size,sizeof(double));
    solver->rhsref  = (double*) calloc (solver->nvars*size,sizeof(double));
    solver->rhs     = (double*) calloc (solver->nvars*size,sizeof(double));
  } else solver->u0 = solver->uref = solver->rhsref = solver->rhs = NULL;
#endif
  solver->hyp     = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->par     = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->source  = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->iblank  = (double*) calloc (size              ,sizeof(double));
  /* grid */
  size = 0;
  for (i=0; i<solver->ndims; i++) size += (solver->dim_local[i]+2*solver->ghosts);
  solver->x     = (double*) calloc (size,sizeof(double));
  solver->dxinv = (double*) calloc (size,sizeof(double));
  /* arrays needed to compute fluxes */
  size = 1;  for (i=0; i<solver->ndims; i++) size *= (solver->dim_local[i]+2*solver->ghosts);
  solver->uC     = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->fluxC  = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->Deriv1 = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->Deriv2 = (double*) calloc (solver->nvars*size,sizeof(double));
  size = 1;  for (i=0; i<solver->ndims; i++) size *= (solver->dim_local[i]+1);
  solver->fluxI = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->uL    = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->uR    = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->fL    = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->fR    = (double*) calloc (solver->nvars*size,sizeof(double));
  /* allocate MPI send/receive buffer arrays */
  int bufdim[solver->ndims], maxbuf = 0;
  for (d = 0; d < solver->ndims; d++) {
    bufdim[d] = 1;
    for (i = 0; i < solver->ndims; i++) {
      if (i == d) bufdim[d] *= solver->ghosts;
      else        bufdim[d] *= solver->dim_local[i];
    }
    if (bufdim[d] > maxbuf) maxbuf = bufdim[d];
  }
  maxbuf *= solver->nvars;
  mpi->maxbuf  = maxbuf;
  mpi->sendbuf = (double*) calloc (2*solver->ndims*maxbuf,sizeof(double));
  mpi->recvbuf = (double*) calloc (2*solver->ndims*maxbuf,sizeof(double));
  /* allocate the volume and boundary integral arrays */
  solver->VolumeIntegral        = (double*) calloc (solver->nvars  ,sizeof(double));
  solver->VolumeIntegralInitial = (double*) calloc (solver->nvars  ,sizeof(double));
  solver->StageBoundaryIntegral = (double*) calloc (2*solver->ndims*solver->nvars,sizeof(double));
  solver->StepBoundaryIntegral  = (double*) calloc (2*solver->ndims*solver->nvars,sizeof(double));
  solver->TotalBoundaryIntegral = (double*) calloc (solver->nvars,sizeof(double));
  solver->ConservationError     = (double*) calloc (solver->nvars,sizeof(double));

  /* initialize function call counts to zero */
  solver->count_hyp = solver->count_par = solver->count_sou = 0;
#ifdef with_petsc
  solver->count_RHSFunction = solver->count_IFunction
    = solver->count_IJacobian = solver->count_IJacFunction 
    = 0;
#endif

  /* Initialize iblank to 1*/
  _ArraySetValue_(solver->iblank,solver->npoints_local_wghosts,1);

  return(0);
}
