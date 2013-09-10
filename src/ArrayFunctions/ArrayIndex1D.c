#include <arrayfunctions.h>

inline int ArrayIndex1D(int N,int *imax,int *i,int *offset,int ghost)
{
  int n,index = i[0];
  for (n = 1; n < N; n++) index = ((index*(imax[n]+2*ghost)) + (i[n]+(offset?offset[n]:0)));
  return(index);
}
