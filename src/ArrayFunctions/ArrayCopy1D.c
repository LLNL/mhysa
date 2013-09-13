#include <arrayfunctions.h>

inline int ArrayCopy1D_double(double *x,double *y,int size)
{
  if (!x || !y) return(1);
  int i;
  for (i = 0; i < size; i++) y[i] = x[i];
  return(0);
}

inline int ArrayCopy1D_int(int *x,int *y,int size)
{
  if (!x || !y) return(1);
  int i;
  for (i = 0; i < size; i++) y[i] = x[i];
  return(0);
}
