#include <arrayfunctions.h>

inline int ArraySetValue_double(double *x,int size,double value)
{
  int n;
  for (n = 0; n < size; n++)  x[n] = value;
  return(0);
}

inline int ArraySetValue_int(int *x,int size,int value)
{
  int n;
  for (n = 0; n < size; n++)  x[n] = value;
  return(0);
}
