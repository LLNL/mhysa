/*! @file MPIBroadcast.c
    @brief Functions to broadcast over all MPI ranks 
    @author Debojyoti Ghosh
*/

#include <mpivars.h>

/*! Broadcast an array of type \a double to all MPI ranks */
int MPIBroadcast_double(
                          double  *x,     /*!< array to broadcast to all ranks */
                          int     size,   /*!< size of array to broadcast */
                          int     root,   /*!< rank from which to broadcast */
                          void    *comm   /*!< MPI communicator within which to broadcast */
                       )
{
#ifndef serial
  MPI_Bcast(x,size,MPI_DOUBLE,root,*((MPI_Comm*)comm));
#endif
  return(0);
}

/*! Broadcast an array of type \a int to all MPI ranks */
int MPIBroadcast_integer(
                          int   *x,     /*!< array to broadcast to all ranks */
                          int   size,   /*!< size of array to broadcast */
                          int   root,   /*!< rank from which to broadcast */
                          void  *comm   /*!< MPI communicator within which to broadcast */
                        )
{
#ifndef serial
  MPI_Bcast(x,size,MPI_INT,root,*((MPI_Comm*)comm));
#endif
  return(0);
}

/*! Broadcast an array of type \a char to all MPI ranks */
int MPIBroadcast_character(
                            char  *x,   /*!< array to broadcast to all ranks */
                            int   size, /*!< size of array to broadcast */
                            int   root, /*!< rank from which to broadcast */
                            void  *comm /*!< MPI communicator within which to broadcast */
                          )
{
#ifndef serial
  MPI_Bcast(x,size,MPI_CHAR,root,*((MPI_Comm*)comm));
#endif
  return(0);
}
