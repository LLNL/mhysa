#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>

inline double ArraySumSquarenD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; ArraySetValue_int(index,ndims,0);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    int v; for (v=0; v<nvars; v++) sum += (x[p*nvars+v]*x[p*nvars+v]);
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(sum);
}
