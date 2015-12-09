/*! @file MPIMin.c
    @brief Functions to compute the minimum over all MPI ranks
    @author Debojyoti Ghosh
*/

#include <mpivars.h>

/*!
  Compute the global minimum over all MPI ranks in a given communicator for
  \a int datatype. 
  + If \a var is an array of size greater than 1, \a global will be an array
    of the same size with each element as the minimum value of he that element 
    in \a var on all the MPI ranks in the given communicator.
*/
int MPIMin_integer(
                    int   *global, /*!< array to contain the global minimums */
                    int   *var,    /*!< the local array */
                    int   size,    /*!< size of the local array */
                    void  *comm    /*!< MPI communicator */
                  )
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_INT,MPI_MIN,*((MPI_Comm*)comm));
#endif
  return(0);
}

/*!
  Compute the global minimum over all MPI ranks in a given communicator for
  \a double datatype. 
  + If \a var is an array of size greater than 1, \a global will be an array
    of the same size with each element as the minimum value of that element 
    in \a var on all the MPI ranks in the given communicator.
*/
int MPIMin_double(
                    double  *global, /*!< array to contain the global minimums */
                    double  *var,    /*!< the local array */ 
                    int     size,    /*!< size of the local array */
                    void    *comm    /*!< MPI communicator */
                 )
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_DOUBLE,MPI_MIN,*((MPI_Comm*)comm));
#endif
  return(0);
}
