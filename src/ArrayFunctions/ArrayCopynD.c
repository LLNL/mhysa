#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>

inline int ArrayCopynD(int ndims,double *x,double *y,int *dim,int g1,int g2,int *index,int nvars)
{
  int ierr = 0;
  if (!y) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"y\" not allocated.\n");
    return(1);
  }
  if (!x) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"x\" not allocated.\n");
    return(1);
  }
  int done = 0;
  ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int p1 = ArrayIndex1D(ndims,dim,index,NULL,g1); 
    int p2 = ArrayIndex1D(ndims,dim,index,NULL,g2);
    int n;
    for (n=0; n<nvars; n++) y[p2*nvars+n] = x[p1*nvars+n];
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(0);
}
