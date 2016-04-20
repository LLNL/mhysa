/*! @file MPIIOGroups.c
    @brief Create I/O groups of MPI ranks
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpivars.h>

/*!
  Create I/O groups of MPI ranks: A scalable approach to file I/O when
  running simulations on a large number of processors (>10,000) is
  partitioning all the MPI ranks into I/O group. Each group has a "leader"
  that:
  + Input - reads the local data of each member rank from the file and sends
    it to that member.
  + Output - gets the local data of each member rank and writes it to a file.

  The number of I/O groups (and hence, the number of I/O ranks reading and 
  writing to files) is specified through #MPIVariables::N_IORanks. Ideally,
  this would correspond to the number of I/O nodes available for the 
  total number of compute nodes being used on a HPC platform. 

  Two extreme cases are:
  + Number of I/O ranks is 1, i.e., just one rank (typically rank 0) is
    responsible for reading and writing the local data of every rank.
  + Number of I/O ranks is equal to the total number of MPI ranks, i.e.,
    each MPI rank reads and writes from and to its own files.

  Neither of the extreme cases are scalable.

  Notes:
  + If the total number of MPI ranks (#MPIVariables::nproc) is not an integer
    multiple of the specified number of I/O groups (#MPIVariables::N_IORanks),
    then only 1 rank is used.
*/
int MPICreateIOGroups(void *m /*!< MPI object of type #MPIVariables*/)
{
  MPIVariables *mpi = (MPIVariables*) m;

#ifndef serial

  int nproc         = mpi->nproc;
  int rank          = mpi->rank;
  int N_IORanks     = mpi->N_IORanks;

  int GroupSize;
  if (nproc%N_IORanks==0) GroupSize = nproc/N_IORanks;
  else {
    if (!rank) {
      printf("Note: in MPICreateIOGroups() - Number of ");
      printf("ranks (nproc) is not divisible by number ");
      printf("of IO ranks (N_IORanks). Readjusting number ");
      printf("of IO ranks to fix this.\n");
    }
    N_IORanks = 1;
    if (!rank) {
      printf("Number of IO Ranks: %d\n",N_IORanks);
    }
    GroupSize = nproc;
  }

  mpi->CommGroup  = rank/GroupSize;
  mpi->IORank     = mpi->CommGroup * GroupSize;

  /* set flag for whether this rank does file I/O */
  if (rank == mpi->IORank)  mpi->IOParticipant = 1;
  else                      mpi->IOParticipant = 0;

  /* save the first and last process of this group */
  mpi->GroupStartRank = mpi->IORank;
  mpi->GroupEndRank   = (mpi->CommGroup+1)*GroupSize;

  /* create a new communicator with the IO participants */
  int i,*FileIORanks;
  MPI_Group WorldGroup, IOGroup;
  FileIORanks = (int*) calloc (N_IORanks,sizeof(int));
  for (i=0; i<N_IORanks; i++) FileIORanks[i] = i*GroupSize;
  MPI_Comm_group(mpi->world,&WorldGroup);
  MPI_Group_incl(WorldGroup,N_IORanks,FileIORanks,&IOGroup);
  MPI_Comm_create(mpi->world,IOGroup,&mpi->IOWorld);
  MPI_Group_free(&IOGroup);
  MPI_Group_free(&WorldGroup);
  free(FileIORanks);

#endif

  return(0);
}
