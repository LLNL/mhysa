/*! @file MPICommunicators.c
    @brief Functions to create and destroy MPI subcommunicators
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Create subcommunicators from MPI_WORLD, where each subcommunicator contains
  MPI ranks along a spatial dimension. Consider a two-dimensional problem, 
  partitioned on 21 MPI ranks as follows:
  @image html mpi_ranks.png
  @image latex mpi_ranks.eps width=0.9\textwidth

  This function will create 10 subcommunicators with the following ranks:
  + 0,1,2,3,4,5,6
  + 7,8,9,10,11,12,13
  + 14,15,16,17,18,19,20
  + 0,7,14
  + 1,8,15
  + 2,9,16
  + 3,10,17
  + 4,11,18
  + 5,12,19
  + 6,13,20

  These subcommunicators are useful for parallel computations along 
  grid lines. For example, a compact finite-difference scheme solves
  implicit systems along grid lines in every spatial dimension. Thus,
  the subcommunicator may be passed on to the parallel systems solver
  instead of MPI_WORLD.
*/
int MPICreateCommunicators(
                            int   ndims,  /*!< Number of spatial dimensions */
                            void  *m      /*!< MPI object of type #MPIVariables */
                          )
{
  MPIVariables *mpi = (MPIVariables*) m;
#ifdef serial
  mpi->comm = NULL;
#else
  int          i,n,color,key;
  int          *ip,*iproc;

  mpi->comm = (MPI_Comm*) calloc (ndims, sizeof(MPI_Comm));
  if (ndims == 1) MPI_Comm_dup(mpi->world,mpi->comm);
  else {
    ip    = (int*) calloc (ndims-1,sizeof(int));
    iproc = (int*) calloc (ndims-1,sizeof(int));
    for (n=0; n<ndims; n++) {
      int tick=0; 
      for (i=0; i<ndims; i++) {
        if (i != n) {
          ip[tick]    = mpi->ip[i];
          iproc[tick] = mpi->iproc[i];
          tick++;
        }
      }
      _ArrayIndex1D_(ndims-1,iproc,ip,0,color); 
      key   = mpi->ip[n];
      MPI_Comm_split(mpi->world,color,key,&mpi->comm[n]);
    }
    free(ip);
    free(iproc);
  }
#endif
  return(0);
}

/*!
  Free the subcommunicators created in MPICreateCommunicators().
*/
int MPIFreeCommunicators(
                          int   ndims,  /*!< Number of spatial dimensions */
                          void  *m      /*!< MPI object of type #MPIVariables */
                        )
{
#ifndef serial
  MPIVariables *mpi = (MPIVariables*) m;
  int          n;
  for (n=0; n<ndims; n++) MPI_Comm_free(&mpi->comm[n]);
  free(mpi->comm);
  if (mpi->IOParticipant) MPI_Comm_free(&mpi->IOWorld);
#endif
  return(0);
}
