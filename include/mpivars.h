#ifndef serial
#include <mpi.h>
#endif

typedef struct mpi_variables {
  int   rank;     /* process rank                                       */
  int   nproc;    /* total number of processes                          */
  int   *iproc;   /* number of processes along each dimension           */
  int   *ip;      /* process rank along each dimension                  */
  int   *is,*ie;  /* global start and end indices along each dimension  */
  int   *bcperiodic; /* flag for periodic BCs along any dimension       */
#ifdef serial
  int   world;
  int   *comm;
#else
  MPI_Comm  world;   /* communicator for all processes                  */
  MPI_Comm  *comm;   /* sub-communicators                               */
#endif
} MPIVariables;

/* Functions */
int MPIBroadcast_double     (double*,int,int,void*);
int MPIBroadcast_integer    (int*,int,int,void*);
int MPIBroadcast_character  (char*,int,int,void*);
int MPICreateCommunicators  (int,void*);
int MPIExchangeBoundaries1D (void*,double*,int,int,int,int);
int MPIExchangeBoundariesnD (int,int,int*,int,void*,double*);
int MPIFreeCommunicators    (int,void*);
int MPIGatherArray1D        (void*,double*,double*,int,int,int,int); 
int MPIGatherArraynD        (int,void*,double*,double*,int*,int*,int,int);
int MPIGetArrayDatanD       (double*,double*,int*,int*,int*,int*,int,int,int,void*);
int MPILocalDomainLimits    (int,int,void*,int*,int*,int*);
int MPIMax_integer          (int*,int*,int,void*);
int MPIMax_double           (double*,double*,int,void*);
int MPIPartition1D          (int,int,int);
int MPIPartitionArraynD     (int,void*,double*,double*,int*,int*,int,int); 
int MPIPartitionArray1D     (void*,double*,double*,int,int,int,int); 
int MPIRank1D               (int,int*,int*);
int MPIRanknD               (int,int,int*,int*);
int MPISum_double           (double*,double*,int,void*);
int MPISum_integer          (int*,int*,int,void*);
