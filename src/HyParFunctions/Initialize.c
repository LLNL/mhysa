#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

int Initialize(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i;

  /* allocations */
  mpi->ip           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->is           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->ie           = (int*) calloc (solver->ndims,sizeof(int));
  solver->dim_local = (int*) calloc (solver->ndims,sizeof(int));

#ifndef serial
  int ierr = 0;

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
  ierr = MPIRanknD(solver->ndims,mpi->rank,mpi->iproc,mpi->ip); CHECKERR(ierr);

  /* calculate local domain sizes along each dimension */
  for (i=0; i<solver->ndims; i++) 
    solver->dim_local[i] = MPIPartition1D(solver->dim_global[i],mpi->iproc[i],mpi->ip[i]);

  /* calculate local domain limits in terms of global domain */
  ierr = MPILocalDomainLimits(solver->ndims,mpi->rank,mpi,solver->dim_global,mpi->is,mpi->ie);
  CHECKERR(ierr);

  /* create sub-communicators for parallel computations along grid lines in each dimension */
  ierr = MPICreateCommunicators(solver->ndims,mpi); CHECKERR(ierr);

#else

  for (i=0; i<solver->ndims; i++) {
    mpi->ip[i]            = 0;
    solver->dim_local[i]  = solver->dim_global[i];
    mpi->iproc[i]         = 1;
    mpi->is[i]            = 0;
    mpi->ie[i]            = solver->dim_local[i];
  }

#endif

  solver->npoints_global = solver->npoints_local = 1;
  for (i=0; i<solver->ndims; i++) solver->npoints_global *= solver->dim_global[i];
  for (i=0; i<solver->ndims; i++) solver->npoints_local  *= solver->dim_local [i];

  /* Allocations */
  if (!mpi->rank) printf("Allocating data arrays.\n");
  solver->index = (int*) calloc (solver->ndims,sizeof(int));
  int size;
  /* state variable */
  size = 1;
  for (i=0; i<solver->ndims; i++) size *= (solver->dim_local[i]+2*solver->ghosts);
  solver->u       = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->hyp     = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->par     = (double*) calloc (solver->nvars*size,sizeof(double));
  solver->source  = (double*) calloc (solver->nvars*size,sizeof(double));
  /* grid */
  size = 0;
  for (i=0; i<solver->ndims; i++) size += (solver->dim_local[i]+2*solver->ghosts);
  solver->x     = (double*) calloc (size,sizeof(double));
  solver->dxinv = (double*) calloc (size,sizeof(double));

  return(0);
}
