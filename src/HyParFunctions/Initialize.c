#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

int Initialize(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0,i;

  /* communicating basic information to all processes */
  ierr = MPIBroadcast_integer(&solver->ndims,1,0); CHECKERR(ierr);
  ierr = MPIBroadcast_integer(&solver->nvars,1,0); CHECKERR(ierr);
  if (mpi->rank) {
    solver->dim_global = (int*) calloc (solver->ndims,sizeof(int));
    mpi->iproc         = (int*) calloc (solver->ndims,sizeof(int));
  }
  ierr = MPIBroadcast_integer(&solver->dim_global[0],solver->ndims,0); CHECKERR(ierr);
  ierr = MPIBroadcast_integer(&mpi->iproc[0]        ,solver->ndims,0); CHECKERR(ierr);

  /* allocations */
  mpi->ip           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->is           = (int*) calloc (solver->ndims,sizeof(int));
  mpi->ie           = (int*) calloc (solver->ndims,sizeof(int));
  solver->dim_local = (int*) calloc (solver->ndims,sizeof(int));

#ifndef serial

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

  int buffer_size = 7;
  int *buffer;
  buffer = (int*) calloc(buffer_size,sizeof(int));
  if (!mpi->rank) {
    buffer[0] = solver->ghosts;
    buffer[1] = solver->n_iter;
    buffer[2] = solver->hyp_space_scheme;
    buffer[3] = solver->par_space_scheme;
    buffer[4] = solver->time_scheme;
    buffer[5] = solver->screen_op_iter;
    buffer[6] = solver->file_op_iter;
  }
  ierr = MPIBroadcast_integer(buffer,buffer_size,0); CHECKERR(ierr);
  solver->ghosts            = buffer[0];
  solver->n_iter            = buffer[1];
  solver->hyp_space_scheme  = buffer[2];
  solver->par_space_scheme  = buffer[3];
  solver->time_scheme       = buffer[4];
  solver->screen_op_iter    = buffer[5];
  solver->file_op_iter      = buffer[6];
  free(buffer);

  ierr = MPIBroadcast_double(&solver->dt,1,0);                               CHECKERR(ierr);
  ierr = MPIBroadcast_character(solver->op_file_format,_MAX_STRING_SIZE_,0); CHECKERR(ierr);
  ierr = MPIBroadcast_character(solver->op_overwrite  ,_MAX_STRING_SIZE_,0); CHECKERR(ierr);

  /* calculate ndims-D rank of each process (ip[]) from rank in MPI_COMM_WORLD */
  ierr = MPIRanknD(solver->ndims,mpi->rank,mpi->iproc,mpi->ip); CHECKERR(ierr);

  /* calculate local domain sizes along each dimension */
  for (i=0; i<solver->ndims; i++) 
    solver->dim_local[i] = MPIPartition1D(solver->dim_global[i],mpi->iproc[i],mpi->ip[i]);

  /* calculate local domain limits in terms of global domain */
  ierr = MPILocalDomainLimits(solver->ndims,mpi->rank,mpi,solver->dim_global,mpi->is,mpi->ie);
  CHECKERR(ierr);

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
  solver->u = (double*) calloc (solver->nvars*size,sizeof(double));
  /* grid */
  size = 0;
  for (i=0; i<solver->ndims; i++) size += solver->dim_local[i];
  solver->x = (double*) calloc (size,sizeof(double));

  return(0);
}
