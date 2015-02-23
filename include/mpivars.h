/*! @file mpivars.h
    @brief MPI related structure and function definitions.
    @author Debojyoti Ghosh
 */

#ifndef serial
#include <mpi.h>
#endif

/*! \def MPIVariables
 * Structure of MPI-related variables.
*/

/*! Structure of MPI-related variables. */
typedef struct mpi_variables {
  int   rank;     /*!< process rank                                       */
  int   nproc;    /*!< total number of processes                          */
  int   *iproc;   /*!< number of processes along each dimension           */
  int   *ip;      /*!< process rank along each dimension                  */
  int   *is,      /*!< global start index of local domain along each dimension  */
        *ie;      /*!< global end index of local domain along each dimension  */
  int   *bcperiodic; /*!< flag for periodic BCs along any dimension       */

#ifdef serial
  int   world; /*!< dummy variable */
  int   *comm; /*!< dummy variable */
#else
  MPI_Comm  world;   /*!< communicator for all processes                  */
  MPI_Comm  *comm;   /*!< sub-communicators                               */
#endif

  int N_IORanks;      /*!< number of IO ranks                             */
  int IOParticipant;  /*!< whether this rank will handle file I/O         */
  int CommGroup;      /*!< I/O group this rank is a part of               */
  int IORank       ;  /*!< Rank of the process this rank will get I/O from*/
  int GroupStartRank; /*!< Starting rank of the IO group                  */
  int GroupEndRank;   /*!< Last rank of the IO group                      */
  MPI_Comm IOWorld;   /*!< Communicator of processes participating in file I/O */

  double *sendbuf, /*!< buffer to send data */
         *recvbuf; /*!< buffer to receive data */
  int    maxbuf;   /*!< maximum buffer size */

} MPIVariables;

/*! broadcast a double to all ranks */
int MPIBroadcast_double     (double*,int,int,void*);
/*! broadcast an integer to all ranks */
int MPIBroadcast_integer    (int*,int,int,void*);
/*! broadcast a character to all ranks */
int MPIBroadcast_character  (char*,int,int,void*);

/*! create communicators required for the tridiagonal solver in compact schemes */
int MPICreateCommunicators  (int,void*);
/*! free/destroy communicators created */
int MPIFreeCommunicators    (int,void*);

/*! create I/O groups for file reading and writing -- group leaders read and writ
 * from/to files and send/receive the data to other ranks in the group */
int MPICreateIOGroups       (void*);

/*! exchange boundary (ghost point) values for an essentially 1D array (like grid
 * coordinates */
int MPIExchangeBoundaries1D (void*,double*,int,int,int,int);

/*! exchange boundary (ghost point) values for an n-dimensional array (like the 
 * solution array) */
int MPIExchangeBoundariesnD (int,int,int*,int,void*,double*);

/*! gather local arrays into a global array for an essentially 1D array */
int MPIGatherArray1D        (void*,double*,double*,int,int,int,int); 
/*! gather local arrays into a global array for an n-dimensional array */
int MPIGatherArraynD        (int,void*,double*,double*,int*,int*,int,int);
/*! partition a global array into local arrays for an n-dimensional array */
int MPIPartitionArraynD     (int,void*,double*,double*,int*,int*,int,int); 
/*! partition a global array into local arrays for an essentially 1D array */
int MPIPartitionArray1D     (void*,double*,double*,int,int,int,int); 

/*! fetch data from an n-dimensional local array on another rank */
int MPIGetArrayDatanD       (double*,double*,int*,int*,int*,int*,int,int,int,void*);

/*! calculate the local domain limits/extend in terms of the global domain */
int MPILocalDomainLimits    (int,int,void*,int*,int*,int*);

/*! find the maximum in an integer array over all ranks */
int MPIMax_integer          (int*,int*,int,void*);
/*! find the maximum in a double array over all ranks */
int MPIMax_double           (double*,double*,int,void*);
/*! find the minimum in an integer array over all ranks */
int MPIMin_integer          (int*,int*,int,void*);
/*! find the minimum in a double array over all ranks */
int MPIMin_double           (double*,double*,int,void*);
/*! calculate the sum of an array of doubles over all ranks */
int MPISum_double           (double*,double*,int,void*);
/*! calculate the sum of an array of integers over all ranks */
int MPISum_integer          (int*,int*,int,void*);

/*! partition (along a dimension) the domain given global size and number of ranks */
int MPIPartition1D          (int,int,int);

/*! calculate 1D rank from the n-dimensional rank */
int MPIRank1D               (int,int*,int*);
/*! calculate the n-dimensional rank from the 1D rank */
int MPIRanknD               (int,int,int*,int*);

/*! generate a unique filename given the rank of the process to let that process
 * write to its own file */
void MPIGetFilename         (char*,void*,char*);
