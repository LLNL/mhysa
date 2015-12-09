/*! @file MPIPartition1D.c
    @brief Compute the local size
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <mpivars.h>

/*!
  Given a 1D array of a given global size \a nglobal, and the total number
  of MPI ranks \a nproc on which it will be partitioned, this function 
  computes the size of the local part of the 1D array on \a rank.
*/
int MPIPartition1D(
                    int nglobal,  /*!< Global size */
                    int nproc,    /*!< Total number of ranks */
                    int rank      /*!< Rank */
                  )
{
  int nlocal;
  if (nglobal%nproc == 0) nlocal = nglobal/nproc;
  else {
    if (rank == nproc-1)  nlocal = nglobal/nproc + nglobal%nproc;
    else                  nlocal = nglobal/nproc;
  }
  return(nlocal);
}
