#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>

inline double ArrayMaxnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; ArraySetValue_int(index,ndims,0);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    int v; 
    for (v=0; v<nvars; v++) {
      if (absolute(x[p*nvars+v]) > sum) sum = absolute(x[p*nvars+v]);
    }
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(sum);
}
