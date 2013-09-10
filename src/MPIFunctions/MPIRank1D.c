#include <mpivars.h>

int MPIRank1D(int ndims,int *iproc,int *ip)
{
  int i,rank = 0, term = 1;
  for (i=0; i<ndims; i++) {
    rank += (ip[i]*term);
    term *= iproc[i];
  }
  
  return(rank);
}
