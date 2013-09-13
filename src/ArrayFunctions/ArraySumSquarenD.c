#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>

inline int ArraySumSquarenD(int ndims,int *dim,int ghosts,int *index,double* x)
{
  if (!x || !index || !dim) return(1);
  int    ierr = 0;
  double sum  = 0;
  int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    sum += (x[p]*x[p]);
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(0);
}
