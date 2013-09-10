#ifndef serial
#include <mpi.h>
#endif

typedef struct mpi_variables {
  int rank;     /* process rank                                       */
  int nproc;    /* total number of processes                          */
  int *iproc;   /* number of processes along each dimension           */
  int *ip;      /* process rank along each dimension                  */
  int *is,*ie;  /* global start and end indices along each dimension  */
} MPIVariables;

/* Functions */
int MPIRank1D             (int,int*,int*);
int MPIRanknD             (int,int,int*,int*);
int MPIPartition1D        (int,int*,int*,int*,int*);
int MPIMax_integer        (int*,int*,int);
int MPIMax_double         (double*,double*,int);
int MPIBroadcast_double   (double*,int,int);
int MPIBroadcast_integer  (int*,int,int);
int MPIBroadcast_character(char*,int,int);
