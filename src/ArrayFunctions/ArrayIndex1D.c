#include <arrayfunctions.h>

inline int ArrayIndex1D(int N,int *imax,int *i,int *offset,int ghost)
{
  int n,index = i[N-1]+(offset?offset[N-1]:0);
  for (n = N-2; n > -1; n--) index = ((index*(imax[n]+2*ghost)) + (i[n]+(offset?offset[n]:0)));
  return(index);
}
