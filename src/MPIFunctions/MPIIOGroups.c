#include <stdio.h>
#include <stdlib.h>
#include <mpivars.h>

int MPICreateIOGroups(void *m)
{
  MPIVariables *mpi = (MPIVariables*) m;

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

  return(0);
}
