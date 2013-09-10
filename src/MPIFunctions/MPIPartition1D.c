#include <stdio.h>
#include <mpivars.h>

int MPIPartition1D(int n,int *nglobal,int *nproc,int *rank,int *nlocal)
{
  int i;
  for (i=0; i<n; i++) {
    if (nglobal[i]%nproc[i] == 0) nlocal[i] = nglobal[i]/nproc[i];
    else {
      if (rank[i] == nproc[i]-1)  nlocal[i] = nglobal[i]/nproc[i] + nglobal[i]%nproc[i];
      else                        nlocal[i] = nglobal[i]/nproc[i];
    }
  }
  return(0);
}
