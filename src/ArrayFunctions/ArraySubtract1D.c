#include <arrayfunctions.h>

inline int ArraySubtract1D_int(int *x,int *a,int *b,int size)
{
  if (!x || !a || !b) return(1);
  int i;
  for (i=0; i<size; i++) x[i] = a[i] - b[i];
  return(0);
}

inline int ArraySubtract1D_double(double *x,double *a,double *b,int size)
{
  if (!x || !a || !b) return(1);
  int i;
  for (i=0; i<size; i++) x[i] = a[i] - b[i];
  return(0);
}
