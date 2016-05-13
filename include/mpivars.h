/*! @file mpivars.h
    @brief MPI related structure and function definitions.
    @author Debojyoti Ghosh
 */

#ifndef serial
#include <mpi.h>
#endif

/*! \def MPIVariables
 *  \brief Structure of MPI-related variables.
 * This structure contains all the variables needed for parallel computations 
 * using the MPI library.
*/

/*! \brief Structure of MPI-related variables.
 *
 * This structure contains all the variables needed for parallel computations 
 * using the MPI library.
*/
typedef struct mpi_variables {
  int   rank;     /*!< Process rank                                       */
  int   nproc;    /*!< Total number of processes                          */
  int   *iproc;   /*!< Number of processes along each dimension           */
  int   *ip;      /*!< Process rank along each dimension                  */
  int   *is,      /*!< Global start index of local domain along each dimension  */
        *ie;      /*!< Global end index of local domain along each dimension  */
  int   *bcperiodic; /*!< Flag for periodic BCs along any dimension       */

#ifdef serial
  int   world; /* Dummy variable */
  int   *comm; /* Dummy variable */
#else
  MPI_Comm  world;   /*!< Communicator for all processes                  */
  MPI_Comm  *comm;   /*!< Sub-communicators                               */
#endif

  int N_IORanks;      /*!< Number of IO ranks                             */
  int IOParticipant;  /*!< Whether this rank will handle file I/O         */
  int CommGroup;      /*!< I/O group this rank is a part of               */
  int IORank       ;  /*!< Rank of the process this rank will get I/O from*/
  int GroupStartRank; /*!< Starting rank of the IO group                  */
  int GroupEndRank;   /*!< Last rank of the IO group                      */
#ifndef serial
  MPI_Comm IOWorld;   /*!< Communicator of processes participating in file I/O */
#endif

  double *sendbuf, /*!< Buffer to send data */
         *recvbuf; /*!< Buffer to receive data */
  int    maxbuf;   /*!< Maximum buffer size */

} MPIVariables;

/*! Broadcast a double to all ranks */
int MPIBroadcast_double     (double*,int,int,void*);
/*! Broadcast an integer to all ranks */
int MPIBroadcast_integer    (int*,int,int,void*);
/*! Broadcast a character to all ranks */
int MPIBroadcast_character  (char*,int,int,void*);

/*! Create communicators required for the tridiagonal solver in compact schemes */
int MPICreateCommunicators  (int,void*);
/*! Free/destroy communicators created */
int MPIFreeCommunicators    (int,void*);

/*! Create I/O groups for file reading and writing -- Group leaders gather data from
 * all other ranks in the group and write to file, and read from file and sends the
 * data to the appropriate ranks in the group. Thus, number of ranks participating
 * in file I/O is equal to the number of groups (which is an input), and can be set
 * to the number of I/O nodes available. */
int MPICreateIOGroups       (void*);

/*! Exchange boundary (ghost point) values for an essentially 1D array (like grid
 * coordinates) */
int MPIExchangeBoundaries1D (void*,double*,int,int,int,int);

/*! Exchange boundary (ghost point) values for an n-dimensional array (like the 
 * solution array) */
int MPIExchangeBoundariesnD (int,int,int*,int,void*,double*);

/*! Gather local arrays into a global array for an essentially 1D array */
int MPIGatherArray1D        (void*,double*,double*,int,int,int,int); 
/*! Gather local arrays into a global array for an n-dimensional array */
int MPIGatherArraynD        (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an n-dimensional array */
int MPIPartitionArraynD     (int,void*,double*,double*,int*,int*,int,int); 
/*! Partition a global array into local arrays for an essentially 1D array */
int MPIPartitionArray1D     (void*,double*,double*,int,int,int,int); 

/*! fetch data from an n-dimensional local array on another rank */
int MPIGetArrayDatanD       (double*,double*,int*,int*,int*,int*,int,int,int,void*);

/*! Calculate the local domain limits/extend in terms of the global domain */
int MPILocalDomainLimits    (int,int,void*,int*,int*,int*);

/*! Find the maximum in an integer array over all ranks */
int MPIMax_integer          (int*,int*,int,void*);
/*! Find the maximum in a long integer array over all ranks */
int MPIMax_long             (long*,long*,int,void*);
/*! Find the maximum in a double array over all ranks */
int MPIMax_double           (double*,double*,int,void*);
/*! Find the minimum in an integer array over all ranks */
int MPIMin_integer          (int*,int*,int,void*);
/*! Find the minimum in a double array over all ranks */
int MPIMin_double           (double*,double*,int,void*);
/*! Calculate the sum of an array of doubles over all ranks */
int MPISum_double           (double*,double*,int,void*);
/*! Calculate the sum of an array of integers over all ranks */
int MPISum_integer          (int*,int*,int,void*);

/*! Partition (along a dimension) the domain given global size and number of ranks */
int MPIPartition1D          (int,int,int);

/*! Calculate 1D rank from the n-dimensional rank */
int MPIRank1D               (int,int*,int*);
/*! Calculate the n-dimensional rank from the 1D rank */
int MPIRanknD               (int,int,int*,int*);

/*! Generate a unique filename given the rank of the process to let that process
 * write to its own file */
void MPIGetFilename         (char*,void*,char*);
