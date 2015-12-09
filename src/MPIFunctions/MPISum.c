/*! @file MPISum.c
    @brief Functions to compute the sum across MPI ranks
    @author Debojyoti Ghosh
*/
#ifndef serial
#include <mpi.h>
#endif

/*!
  Compute the global sum over all MPI ranks in a given communicator for
  \a int datatype. 
  + If \a var is an array of size greater than 1, \a global will be an array
    of the same size with each element as the sum of that element 
    in \a var on all the MPI ranks in the given communicator.
*/
int MPISum_integer(
                    int   *global, /*!< array to contain the global sums */
                    int   *var,    /*!< the local array */
                    int   size,    /*!< size of the local array */
                    void  *comm    /*!< MPI communicator */
                  )
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_INT,MPI_SUM,*((MPI_Comm*)comm));
#endif
  return(0);
}

/*!
  Compute the global sum over all MPI ranks in a given communicator for
  \a double datatype. 
  + If \a var is an array of size greater than 1, \a global will be an array
    of the same size with each element as the sum of that element 
    in \a var on all the MPI ranks in the given communicator.
*/
int MPISum_double(
                    double  *global, /*!< array to contain the global sums */
                    double  *var,    /*!< the local array */
                    int     size,    /*!< size of the local array */
                    void    *comm    /*!< MPI communicator */
                 )
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_DOUBLE,MPI_SUM,*((MPI_Comm*)comm));
#endif
  return(0);
}
