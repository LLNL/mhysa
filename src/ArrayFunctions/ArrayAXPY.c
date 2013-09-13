#include <arrayfunctions.h>

inline int ArrayAXPY(double *x,double a,double *y,int size)
{
  if (!x || !y) return(1);
  int i;
  for (i=0; i<size; i++) y[i] += a*x[i];
  return(0);
}

