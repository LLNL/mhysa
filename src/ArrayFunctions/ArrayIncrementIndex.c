#include <arrayfunctions.h>

inline int ArrayIncrementIndex(int N,int *imax,int *i)
{
  int n = 0;
  while (n < N) {
    if (i[n] == imax[n]-1) {
      i[n] = 0;
      n++;
    } else {
      i[n]++;
      break;
    }
  }
  if (n == N) return(1);
  else        return(0);
}
